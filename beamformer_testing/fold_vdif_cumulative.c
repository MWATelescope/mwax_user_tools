// Fold a file of complex samples at a specified period and DM
// - mainly aimed at VDIF files, but could be raw data or CODIF too
// - channelises, detects power of both polarisations and sums them before folding
// - assumes total number of samples to fold is 8 seconds (one MWA subfile)
// - assumes the two pols are interleaved
// - only 8 bit supported at present
// - saves a cumulative profile that is updated each time the program is run, summing
//   the profile for the current subfile with a loaded profile from past runs
// - the MJD seconds from the first frame of the input VDIF file is used to
//   determine the folding phase offset, so input VDIF files do not have to be contiguous
// - to obtain the profile for just the current subfile, specify an all-zeros input profile
//
// Build: gcc -o fold_vdif_cumulative fold_vdif_cumulative.c -lfftw3 -lm -lvdifio
//
// Ian Morrison December 2025

// system include files
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <vdifio.h>

#define TOTAL_SAMPLES 10240000                // total number of complex samples in each pol in an 1.28 MHz 8 second MWA subfile
#define FFT_SIZE 200                          // FFT length for channelisation
#define SAMPLE_RATE 1280000.0                 // sample rate in Hz for 1.28 MHz bandwidth

#define DEFAULT_HEADER_SIZE 32                // default is VDIF
#define DEFAULT_DATA_FRAME_SIZE 8192          // 8192 bytes
#define DEFAULT_NUM_INPUT_BITS 8              // twice this for a complex sample
#define DEFAULT_CENTRE_FREQUENCY_MHZ 150.0    // centre frequency in MHz
#define DEFAULT_MJD_START_SECONDS 5200000000L //
#define DEFAULT_NUM_PHASE_BINS 256            // number of phase bins for folding the data
#define DEFAULT_PROFILE_FILENAME "pulse_profile.dat"

void usage()
{
  fprintf(stdout, "\fold_vdif - fold a VDIF file at a specified period and DM\n\n"
    "  Usage: fold_vdif <options>\n"
    "         -H <int>          header size in bytes (0 for raw data, 32 for VDIF, 64 for CODIF) [default %d]\n"
    "         -d <int>          data frame size in bytes [default %d]\n"
    "         -b <int>          bits per sample (real and imag) [default %d]\n"
    "         -u                if this option is present, data is assumed to be unsigned (default is signed)\n"
    "         -B <int>          number of phase bins for folding the data [default %d]\n"
    "         -P <period>       period at which to fold the data (diagnostic) [default no folding]\n"
    "         -D <num>          DM for incoherent de-dispersion [default no de-dispersion]\n"
    "         -M <num>          MJD seconds for the start of folding (must be earlier than the MJS seconds of any VDIF file input) [default %lu]\n"
    "         -f <freq>         centre frequency in MHz (for de-dispersion) [default %f MHz]\n"
    "         -l                output the pulse profile in log scale (dB)\n"
    "         -n <int>          number of frames to process - minimum 1 [default all frames in the input file]\n"
    "         -N                normalise the output pulse profile peak to 1.0 (0.0 for log output)\n"
    "         -i <filename>     input VDIF filename"
    "         -p <filename>     input/output pulse profile filename [default %s]\n"
    "         -h                print this usage information\n\n",
    DEFAULT_HEADER_SIZE, DEFAULT_DATA_FRAME_SIZE, DEFAULT_NUM_INPUT_BITS, DEFAULT_NUM_PHASE_BINS, DEFAULT_MJD_START_SECONDS, DEFAULT_CENTRE_FREQUENCY_MHZ, DEFAULT_PROFILE_FILENAME);
}


int main(int argc, char *argv[])
{
  FILE *fp_vdif;
  FILE *fp_profile;

  // command line option parameters
  struct option options[] = {
    {"H",           required_argument, 0, 'H'},
    {"d",           required_argument, 0, 'd'},
    {"b",           required_argument, 0, 'b'},
    {"u",           no_argument,       0, 'u'},
    {"B",           required_argument, 0, 'B'},
    {"P",           required_argument, 0, 'P'},
    {"D",           required_argument, 0, 'D'},
    {"M",           required_argument, 0, 'M'},
    {"f",           required_argument, 0, 'f'},
    {"l",           no_argument,       0, 'l'},
    {"n",           required_argument, 0, 'n'},
    {"N",           no_argument,       0, 'N'},
    {"i",           required_argument, 0, 'i'},
    {"p",           required_argument, 0, 'o'},
    {"h",           no_argument,       0, 'h'},
    {0, 0, 0, 0}
  };

  int opt, ss;
  char vdif_fname[200];
  vdif_fname[0]= '\0';  // empty string
  char profile_fname[200] = DEFAULT_PROFILE_FILENAME;

  int log_output = 0;  // linear power output by default
  int frames_to_process = -1;  // process all frames by default
  int normalise_output = 0;  // do not normalise the output by default 
  int frame_hsize = DEFAULT_HEADER_SIZE;
  int frame_dsize = DEFAULT_DATA_FRAME_SIZE;
  int bits_per_sample = DEFAULT_NUM_INPUT_BITS;  // 2x this for full complex sample
  int format = 0;              // 0 = signed (default), 1 = unsigned
  int num_phase_bins = DEFAULT_NUM_PHASE_BINS;
  double folding_period = 1.0;  // default folding period to prevent divide by zero
  double folding_dm = 0.0;
  uint64_t folding_start_mjd_seconds = DEFAULT_MJD_START_SECONDS;
  double centre_frequency_MHz = DEFAULT_CENTRE_FREQUENCY_MHZ;

  // parse command line options
  while (1)
  {
    opt = getopt_long_only(argc, argv, "H:d:b:uB:P:D:M:f:ln:Ni:p:h", options, NULL);
    
    if (opt == EOF) break;
    
    switch (opt)
    {
      case 'H':
        frame_hsize = atoi(optarg);
        if ((frame_hsize != 0) && (frame_hsize != 32) && (frame_hsize != 64))
        {
          fprintf(stderr, "ERROR: invalid header size %d (must be 0, 32 or 64)\n", frame_hsize);
          exit(EXIT_FAILURE);
        }
        break;

      case 'd':
        frame_dsize = atoi(optarg);
        break;

      case 'b':
        bits_per_sample = atoi(optarg);
        if (bits_per_sample != 8)
        {
          fprintf(stderr, "ERROR: invalid bits per sample %d (must be 8)\n", bits_per_sample);
          exit(EXIT_FAILURE);
        }
        break;

      case 'u':
        format = 1;  // unsigned
        break;

      case 'B':
        num_phase_bins = atoi(optarg);
        if (num_phase_bins < 1)
        {
          fprintf(stderr, "ERROR: bad -B option %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'P':
        folding_period = atof(optarg);
        if (folding_period <= 0.0)
        {
          fprintf(stderr, "ERROR: bad -P option %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'D':
        folding_dm = atof(optarg);
        break;

      case 'M':
        folding_start_mjd_seconds = atoi(optarg);
        if (folding_start_mjd_seconds < DEFAULT_MJD_START_SECONDS)
        {
          fprintf(stderr, "ERROR: bad -M option %s, folding precision will be degraded\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'f':
        centre_frequency_MHz = atof(optarg);
        break;

      case 'l':
        log_output = 1;
        break;

      case 'n':
        ss = sscanf(optarg, "%d", &frames_to_process);
        if (ss!=1 || frames_to_process < 1)
        {
          fprintf(stderr, "ERROR: bad -n option %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'N':
        normalise_output = 1;
        break;

      case 'i':
        ss = sscanf(optarg, "%s", vdif_fname);
        if (ss!=1)
        {
          fprintf(stderr, "ERROR: bad -i option %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'p':
        ss = sscanf(optarg, "%s", profile_fname);
        if (ss!=1)
        {
          fprintf(stderr, "ERROR: bad -p option %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'h':
        usage();
        exit(EXIT_FAILURE);

      default:
        fprintf(stderr, "ERROR: unknown command line option %c\n", opt);
        exit(EXIT_FAILURE);
    }
  }

  int bits_per_frame = 8 * frame_dsize;        // frame_dsize is in bytes
  if (bits_per_frame % (2*bits_per_sample))
  {
    fprintf(stderr, "ERROR: cannot fit an integer number of complex samples (each 2 x %d bits) into a single frame of size %d bits, exiting\n", bits_per_sample, bits_per_frame);
    exit(EXIT_FAILURE);
  }

  // check for an input VDIF file, and if there is one specified, try to open it
  if (vdif_fname[0] == '\0')
  {
    fprintf(stderr, "ERROR: no input VDIF file specified (with the -i option)\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    fp_vdif = fopen(vdif_fname, "r");
    if (fp_vdif == NULL)
    {
      fprintf(stderr, "ERROR: cannot open input file %s\n", vdif_fname);
      exit(EXIT_FAILURE);
    }
  }
  fp_profile = fopen(profile_fname, "r");
  if (fp_profile == NULL)
  {
    fprintf(stderr, "ERROR: cannot open profile file %s\n", profile_fname);
    fclose(fp_vdif);
    exit(EXIT_FAILURE);
  }
  fprintf(stdout, "INFO: input VDIF file = %s\n", vdif_fname);
  fprintf(stdout, "INFO: input/output profile file = %s\n", profile_fname);

  // read in the existing profile data
  float *phase_bins = (float *)malloc(num_phase_bins * sizeof(float));
  for (int i = 0; i < num_phase_bins; i++)
  {
    float bin_value;
    int ret = fscanf(fp_profile, "%f", &bin_value);
    if (ret == 1)
    {
      phase_bins[i] = bin_value;
    }
    else
    {
      fprintf(stderr, "ERROR: cannot read bin %d from profile file %s\n", i, profile_fname);
      fclose(fp_vdif);
      fclose(fp_profile);
      free(phase_bins);
      exit(EXIT_FAILURE);
    }
  }
  fclose(fp_profile);  // close for now, will reopen later for writing

  // set the FFT size
  int samples_per_frame_per_pol = frame_dsize / (2*(bits_per_sample/4));  // number of complex samples per frame per pol
  fprintf(stdout, "INFO: samples per frame per pol is %d\n", samples_per_frame_per_pol);
  int fft_size = FFT_SIZE;
  int samples_per_channel = TOTAL_SAMPLES / fft_size;
  fprintf(stdout, "INFO: channelising into %d channels, %d samples in each channel (number of FFTs)\n", fft_size, samples_per_channel);

  // allocate host memory buffers
  int total_frames = TOTAL_SAMPLES / samples_per_frame_per_pol;
  char *input_frame = (char *)malloc(frame_dsize);
  fftw_complex *input_pol0 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*TOTAL_SAMPLES);
  fftw_complex *input_pol1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*TOTAL_SAMPLES);
  fftw_complex *output_pol0 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*TOTAL_SAMPLES);
  fftw_complex *output_pol1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*TOTAL_SAMPLES);
  double *input_double;
  double *output_double;
  fftw_plan plan_pol0 = fftw_plan_dft_1d(fft_size, input_pol0, output_pol0, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan plan_pol1 = fftw_plan_dft_1d(fft_size, input_pol1, output_pol1, FFTW_FORWARD, FFTW_MEASURE);
  float *power_both_pols = (float *)malloc(TOTAL_SAMPLES * sizeof(float));

  // mask for signed vs unsigned data (used for 8 and 16 bit cases)
  int data_mask;
  if (format == 0)   // signed data
    data_mask = 0;
  else              // unsigned data
    data_mask = 1 << (bits_per_sample - 1);   // will be 0x80 for 8-bit data, 0x8000 for 16-bit data

  // read the input VDIF file and process all available frames
  size_t bytes_read;
  long int frames_read = 0;
  uint64_t first_mjd_secs_this_file;

  for (int frame = 0; frame < total_frames; frame++)
  {
    bytes_read = fread(input_frame, 1, frame_hsize, fp_vdif);    // read the header bytes
    if (frame == 0)
    {
      first_mjd_secs_this_file = getVDIFFrameMJDSec((vdif_header *)input_frame);
      fprintf(stdout, "MJD seconds of first frame of the current VDIF input file = %lu\n", first_mjd_secs_this_file);
    }

    bytes_read = fread(input_frame, 1, frame_dsize, fp_vdif);    // overwrite with the data

    if (frame == 0)
    {
      // print the first few bytes of the first frame
      fprintf(stdout, "First 16 bytes of first data frame:\n");
      for (int i=0; i<16; i++)
      {
        fprintf(stdout, "    byte[%d] = %d\n", i, (char)(input_frame[i] ^ data_mask));
      }
    }

    if (bytes_read != frame_dsize)
    {
      fprintf(stdout, "INFO: cannot read %d bytes from input file, assuming reached EOF\n", frame_dsize);
    }
    else
    {
      frames_read++;
      // copy the input data to the two pol sample buffers
      for (int i = 0; i < samples_per_frame_per_pol; i++)
      {
        input_double = (double *)(input_pol0 + (frame * samples_per_frame_per_pol));
        input_double[2*i] = (double)(char)(input_frame[4*i + 0] ^ data_mask);
        input_double[2*i+1] = (double)(char)(input_frame[4*i + 1] ^ data_mask);
        input_double = (double *)(input_pol1 + (frame * samples_per_frame_per_pol));
        input_double[2*i] = (double)(char)(input_frame[4*i + 2] ^ data_mask);
        input_double[2*i+1] = (double)(char)(input_frame[4*i + 3] ^ data_mask);
      }
    }
  }
  fprintf(stdout, "INFO: number of frames read = %ld\n", frames_read);

#if 1
  // print the first samples of the the input file
  fprintf(stdout, "First complex voltage samples of the input file:\n");
  for (int i=0; i<10; i++)
  {
    input_double = (double *)input_pol0;
    fprintf(stdout, "    pol0 sample[%d] = %f   %f\n", i, input_double[2*i], input_double[2*i+1]);
    input_double = (double *)input_pol1;
    fprintf(stdout, "    pol1 sample[%d] = %f   %f\n", i, input_double[2*i], input_double[2*i+1]);
  }
#endif

  // FFT channelise each pol
  for (int ifft = 0; ifft < samples_per_channel; ifft++)
  {
    fftw_complex *in_ifft  = input_pol0 +  (ifft * fft_size);
    fftw_complex *out_ifft = output_pol0 + (ifft * fft_size);
    fftw_execute_dft(plan_pol0, in_ifft, out_ifft);
    in_ifft  = input_pol1 +  (ifft * fft_size);
    out_ifft = output_pol1 + (ifft * fft_size);
    fftw_execute_dft(plan_pol1, in_ifft, out_ifft);
  }
  fprintf(stdout, "INFO: %d FFTs completed for each pol\n", samples_per_channel);

#if 1
  // print the first samples of the the channelised input file
  fprintf(stdout, "First complex voltage samples of the channelised input file:\n");
  for (int i=0; i<10; i++)
  {
    output_double = (double *)output_pol0;
    fprintf(stdout, "    pol0 channelised[%d] = %f   %f\n", i, output_double[2*i], output_double[2*i+1]);
    output_double = (double *)output_pol1;
    fprintf(stdout, "    pol1 channelised[%d] = %f   %f\n", i, output_double[2*i], output_double[2*i+1]);
  }
#endif

#if 1
  // compute the power of both polarisations and sum them
  fprintf(stdout, "INFO: calculating the power of both polarisations and summing them\n");
  for (int i=0; i<TOTAL_SAMPLES; i++)
  {
    power_both_pols[i] = (float)(creal(output_pol0[i])*creal(output_pol0[i]) + cimag(output_pol0[i])*cimag(output_pol0[i]) +
                                creal(output_pol1[i])*creal(output_pol1[i]) + cimag(output_pol1[i])*cimag(output_pol1[i]));
  }
#endif

#if 1
  // print the first powers
  fprintf(stdout, "First powers:\n");
  for (int i=0; i<FFT_SIZE; i++)   // one whole FFT size worth
  {
    fprintf(stdout, "    power both pols[%d] = %f\n", i, power_both_pols[i]);
  }
#endif

  // at this point we have Stokes I power samples in power_both_pols[] array, channelised

  fprintf(stdout, "Folding each fine channel\n");
  int *num_values_per_bin = (int *)calloc(num_phase_bins, sizeof(int));           // bin counts all initialised to zero
  double *channel_dm_time_offsets = (double *)malloc(fft_size * sizeof(double));  // in seconds
  double max_offset = -10.0;  // set to a large negative number initially
  double fine_chan_width_MHz = (1.0E-6 * (double)SAMPLE_RATE) / (double)fft_size;
  int channel;   // the channel number in ascending sky frequency order
  for (int c = 0; c < fft_size; c++)
  {
    if (c <= (fft_size / 2))
      channel = c;
    else
      channel = c - fft_size;
    double channel_freq_MHz = centre_frequency_MHz + ((double)channel * fine_chan_width_MHz);
    if (channel_freq_MHz != 0.0f)
      channel_dm_time_offsets[c] = 4.148808E3 * folding_dm / (channel_freq_MHz * channel_freq_MHz);  // in seconds (4.148808E3 is dispersion constant MHz^2 cm^3 pc^-1 s)
    else
      channel_dm_time_offsets[c] = 0.0f;
    if (channel_dm_time_offsets[c] > max_offset) max_offset = channel_dm_time_offsets[c];
  }
  // adjust the offsets so that the largest is zero, and all others are negative
  for (int c = 0; c < fft_size; c++)
  {
    channel_dm_time_offsets[c] -= max_offset;
    if (channel_dm_time_offsets[c] > 0.0)   // just in case of rounding error
      channel_dm_time_offsets[c] = 0.0;
  }
  #if 1  // debug
  for (int c = 0; c < fft_size; c++)
  {
    fprintf(stdout, "channel %4d: DM time offset = %.6f s\n", c, channel_dm_time_offsets[c]);
  }
  #endif

  double sample_period = (double)fft_size / SAMPLE_RATE;          // time between successive samples in each channel

  // fold each channel separately, applying incoherent de-dispersion if requested
  int i, j;
  for (i = 0; i < fft_size; i++)
  {
    double start_time = (double)(first_mjd_secs_this_file - folding_start_mjd_seconds) - channel_dm_time_offsets[i];  // in seconds, should not be negative (if MJD offsets are >= 0.0)

    // apply incoherent de-dispersion
    if (start_time < 0.0)  // guarantee not negative (in case of rounding errors)
      start_time = 0.0;

    // fold each sample in this channel
    // - first sample is at start_time, increment by sample_period on each successive sample
    for (j = 0; j < samples_per_channel; j++)
    {
      double fractional_period = fmod(start_time, folding_period) / folding_period;   // range is [0.0, 1.0)], with the value 1.0 possible but unlikely
      int phase_bin = (int)floor(fractional_period * (double)num_phase_bins);         // range is [0, num_phase_bins)], with the value num_phase_bins possible but unlikely
      if (phase_bin == num_phase_bins) phase_bin = num_phase_bins - 1;                // ensures range is [0, num_phase_bins-1]
      num_values_per_bin[phase_bin]++;
      phase_bins[phase_bin] += power_both_pols[j * fft_size + i];

      // next sample
      start_time += sample_period;
    }
  }

  // reopen the profile file for writing the updated profile
  fp_profile = fopen(profile_fname, "w");

  // write the output profile to file
  for (i = 0; i < num_phase_bins; i++)
  {
    fprintf(fp_profile, "%f\n", phase_bins[i]);
  }

  if (normalise_output)
  {
    // find the maximum bin value
    float max_bin_value = -1.0f;
    for (i = 0; i < num_phase_bins; i++)
    {
      if (phase_bins[i] > max_bin_value) max_bin_value = phase_bins[i];
    }
    if (max_bin_value == 0.0) max_bin_value = 1.0;  // don't normalise in zero case (avoids NaN errors)

    // normalise all bins to the maximum value
    float normalisation_factor = 1.0f / max_bin_value;
    for (i = 0; i < num_phase_bins; i++)
    {
      phase_bins[i] *= normalisation_factor;
    }
  }

  // print all the bin values
  fprintf(stdout, "Pulse profile:\n");
  for (i = 0; i < num_phase_bins; i++)
  {
    if (log_output)
    {
      phase_bins[i] = 10.0f * log10f(phase_bins[i] + 1.0e-20f);  // add a small value to avoid log10(0)
    }
    fprintf(stdout, "  bin[%3d]: %f   (num values = %d)\n", i, phase_bins[i], num_values_per_bin[i]);
  }

  cleanup:

  fclose(fp_vdif);
  fclose(fp_profile);

  // free host memory
  free(input_frame);
  free(power_both_pols);
  free(phase_bins);
  free(num_values_per_bin);
  free(channel_dm_time_offsets);

  // destroy the FFTW stuff
  fftw_destroy_plan(plan_pol0);
  fftw_destroy_plan(plan_pol1);
  fftw_free(input_pol0);
  fftw_free(input_pol1);
  fftw_free(output_pol0);
  fftw_free(output_pol1);

  // final message
  fprintf(stdout, "INFO: all done, exiting\n\n");
  return EXIT_SUCCESS;
}
