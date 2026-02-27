// Fold a file of complex samples at a specified period and DM
// - mainly aimed at VDIF files, but could be raw data or CODIF too
// - channelises, detects power of both polarisations and sums them before folding
// - assumes total number of samples to fold is 8 seconds (one MWA subfile)
// - assumes the two pols are interleaved
// - only 8 bit supported at present
//
// Build: gcc -O3 -o fold_vdif fold_vdif.c plot.c -lfftw3 -lm
//
// Ian Morrison November 2025

// system include files
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "plot.h"

#define FFT_SIZE 200          // FFT length for channelisation
#define SAMPLE_RATE 1280000.0 // sample rate in Hz for 1.28 MHz bandwidth

#define DEFAULT_START_TIME 0                // start time in seconds for folding
#define DEFAULT_END_TIME -1                 // end time in seconds for folding (default is -1 to fold all available data)
#define DEFAULT_HEADER_SIZE 32              // default is VDIF
#define DEFAULT_DATA_FRAME_SIZE 8000        // 8000 bytes
#define DEFAULT_NUM_INPUT_BITS 8            // twice this for a complex sample
#define DEFAULT_CENTRE_FREQUENCY_MHZ 150.0f // centre frequency in MHz
#define DEFAULT_NUM_PHASE_BINS 128          // number of phase bins for folding the data
#define DEFAULT_OUTPUT_FILENAME "pulse_profile.dat"

int parse_vdif_filename(const char *path,
                        int *obsid,
                        int *channel,
                        int *beam)
{
  // Step past the last '/'
  const char *fname = strrchr(path, '/');
  fname = fname ? fname + 1 : path;

  // Parse fields using sscanf
  if (sscanf(fname, "%10d_ch%d_beam%d", obsid, channel, beam) != 3)
  {
    fprintf(stderr, "ERROR: could not parse filename '%s'\n", fname);
    return -1;
  }

  return 0;
}

long get_file_size(const char *filename)
{
  struct stat st;
  if (stat(filename, &st) != 0)
    return -1;
  return st.st_size;
}

void usage()
{
  fprintf(stdout, "\nfold_vdif - fold a VDIF file at a specified period and DM\n\n"
                  "  Usage: fold_vdif <options>\n"
                  "         -s <int>          start time in seconds for folding [default %d s]\n"
                  "         -e <int>          end time in seconds for folding [default %d s]\n"
                  "         -H <int>          header size in bytes (0 for raw data, 32 for VDIF, 64 for CODIF) [default %d]\n"
                  "         -d <int>          data frame size in bytes [default %d]\n"
                  "         -b <int>          bits per sample (real and imag) [default %d]\n"
                  "         -u                if this option is present, data is assumed to be unsigned (default is signed)\n"
                  "         -B <int>          number of phase bins for folding the data [default %d]\n"
                  "         -P <period>       period at which to fold the data (diagnostic) [default no folding]\n"
                  "         -D <num>          DM for incoherent de-dispersion [default no de-dispersion]\n"
                  "         -O <offset>       non-negative phase offset time in seconds (to place the pulse at a desired phase) [default 0.0 s]\n"
                  "         -f <freq>         centre frequency in MHz (for de-dispersion) [default %f MHz]\n"
                  "         -l                output the pulse profile in log scale (dB)\n"
                  "         -n                number of frames to process - minimum 1 [default all frames in the input file]\n"
                  "         -N                normalise the output pulse profile peak to 1.0 (0.0 for log output)\n"
                  "         -i                input filename\n"
                  "         -o                output filename [default %s]\n"
                  "         -h                print this usage information\n\n",
          DEFAULT_START_TIME, DEFAULT_END_TIME, DEFAULT_HEADER_SIZE, DEFAULT_DATA_FRAME_SIZE, DEFAULT_NUM_INPUT_BITS, DEFAULT_NUM_PHASE_BINS, DEFAULT_CENTRE_FREQUENCY_MHZ, DEFAULT_OUTPUT_FILENAME);
}

int main(int argc, char *argv[])
{
  FILE *fp_in;
  FILE *fp_out;

  // command line option parameters
  struct option options[] = {
      {"s", required_argument, 0, 's'},
      {"e", required_argument, 0, 'e'},
      {"H", required_argument, 0, 'H'},
      {"d", required_argument, 0, 'd'},
      {"b", required_argument, 0, 'b'},
      {"u", no_argument, 0, 'u'},
      {"B", required_argument, 0, 'B'},
      {"P", required_argument, 0, 'P'},
      {"D", required_argument, 0, 'D'},
      {"O", required_argument, 0, 'O'},
      {"f", required_argument, 0, 'f'},
      {"l", no_argument, 0, 'l'},
      {"n", required_argument, 0, 'n'},
      {"N", no_argument, 0, 'N'},
      {"i", required_argument, 0, 'i'},
      {"o", required_argument, 0, 'o'},
      {"h", no_argument, 0, 'h'},
      {0, 0, 0, 0}};

  int opt, ss;
  char in_fname[200];
  in_fname[0] = '\0'; // empty string
  char out_fname[200] = DEFAULT_OUTPUT_FILENAME;

  int start_time = DEFAULT_START_TIME;
  int end_time = DEFAULT_END_TIME;
  int log_output = 0;         // linear power output by default
  int frames_to_process = -1; // process all frames by default
  int normalise_output = 0;   // do not normalise the output by default
  int frame_hsize = DEFAULT_HEADER_SIZE;
  int frame_dsize = DEFAULT_DATA_FRAME_SIZE;
  int bits_per_sample = DEFAULT_NUM_INPUT_BITS; // 2x this for full complex sample
  int format = 0;                               // 0 = signed (default), 1 = unsigned
  int num_phase_bins = DEFAULT_NUM_PHASE_BINS;
  double folding_period = 1.0; // default folding period to prevent divide by zero
  float folding_dm = 0.0;
  float phase_offset_time = 0.0f; // in seconds
  float centre_frequency_MHz = DEFAULT_CENTRE_FREQUENCY_MHZ;

  // parse command line options
  while (1)
  {
    opt = getopt_long_only(argc, argv, "s:e:H:d:b:uB:P:D:O:f:ln:Ni:o:h", options, NULL);

    if (opt == EOF)
      break;

    switch (opt)
    {
    case 's':
      start_time = atoi(optarg);
      if (start_time < 0)
      {
        fprintf(stderr, "ERROR: -s (start time) must be >=0 %s\n", optarg);
        exit(EXIT_FAILURE);
      }
      break;

    case 'e':
      end_time = atoi(optarg);
      break;

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
      format = 1; // unsigned
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

    case 'O':
      phase_offset_time = atof(optarg);
      if (phase_offset_time < 0.0)
      {
        fprintf(stderr, "ERROR: bad -O option %s\n", optarg);
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
      if (ss != 1 || frames_to_process < 1)
      {
        fprintf(stderr, "ERROR: bad -n option %s\n", optarg);
        exit(EXIT_FAILURE);
      }
      break;

    case 'N':
      normalise_output = 1;
      break;

    case 'i':
      ss = sscanf(optarg, "%s", in_fname);
      if (ss != 1)
      {
        fprintf(stderr, "ERROR: bad -file option %s\n", optarg);
        exit(EXIT_FAILURE);
      }
      break;

    case 'o':
      ss = sscanf(optarg, "%s", out_fname);
      if (ss != 1)
      {
        fprintf(stderr, "ERROR: bad -file option %s\n", optarg);
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

  int bits_per_frame = 8 * frame_dsize; // frame_dsize is in bytes
  if (bits_per_frame % (2 * bits_per_sample))
  {
    fprintf(stderr, "ERROR: cannot fit an integer number of complex samples (each 2 x %d bits) into a single frame of size %d bits, exiting\n", bits_per_sample, bits_per_frame);
    exit(EXIT_FAILURE);
  }

  // Determine total frames
  long int file_size = get_file_size(in_fname);
  if (file_size < 0)
  {
    fprintf(stderr, "ERROR: cannot determine size of input file %s\n", in_fname);
    exit(EXIT_FAILURE);
  }

  int obsid = 0, rec_channel = 0, beam = 0;

  if (parse_vdif_filename(in_fname, &obsid, &rec_channel, &beam) != 0)
  {
    fprintf(stderr, "ERROR: could not parse input filename %s\n", in_fname);
    exit(EXIT_FAILURE);
  }

  int frame_size = frame_hsize + frame_dsize;
  long total_frames = file_size / frame_size;
  // 4 because each complex sample has 2 pols and each pol has real and imag parts, each of bits_per_sample bits
  long samples_per_frame = bits_per_frame / (bits_per_sample * 4);
  long total_samples = total_frames * samples_per_frame;
  int total_seconds = (int)((double)total_samples / SAMPLE_RATE);

  if (end_time > total_seconds)
  {
    fprintf(stderr, "WARNING: end time %d s is greater than total available data duration %d s, adjusting end time to %d s\n", end_time, total_seconds, total_seconds);
    end_time = total_seconds;
  }
  else if (end_time == -1) // default case to fold all available data
  {
    end_time = total_seconds;
  }

  int selected_seconds = end_time - start_time;
  long selected_samples = (long)(selected_seconds * SAMPLE_RATE);
  int selected_start_frame = (int)((double)start_time * SAMPLE_RATE / samples_per_frame);
  int selected_end_frame = (int)((double)end_time * SAMPLE_RATE / samples_per_frame);

  printf("INFO: OBSID: %d\n", obsid);
  printf("INFO: channel: %d\n", rec_channel);
  printf("INFO: beam number: %d\n\n", beam);
  printf("INFO: input file size is %ld bytes\n", file_size);
  printf("INFO: frame size is %d bytes (header %d + data %d)\n", frame_size, frame_hsize, frame_dsize);
  printf("INFO: total frames in the file is %ld\n", total_frames);
  printf("INFO: samples per frame is %ld\n", samples_per_frame);
  printf("INFO: total samples in the file is %ld\n", total_samples);
  printf("INFO: total seconds in the file is %d\n\n", total_seconds);
  printf("INFO: SELECTED seconds in the file is %d\n", selected_seconds);
  printf("INFO: SELECTED samples in the file is %ld\n", selected_samples);
  printf("INFO: SELECTED start frame is %d\n", selected_start_frame);
  printf("INFO: SELECTED end frame is %d\n\n", selected_end_frame);

  // check for an input file, and if there is one specified, try to open it
  if (in_fname[0] == '\0')
  {
    fprintf(stderr, "ERROR: no input file specified (with the -i option)\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    fp_in = fopen(in_fname, "r");
    if (fp_in == NULL)
    {
      fprintf(stderr, "ERROR: cannot open input file %s\n", in_fname);
      exit(EXIT_FAILURE);
    }
  }
  fp_out = fopen(out_fname, "w");
  if (fp_out == NULL)
  {
    fprintf(stderr, "ERROR: cannot open output file %s\n", out_fname);
    fclose(fp_in);
    exit(EXIT_FAILURE);
  }
  fprintf(stdout, "INFO: input file = %s\n", in_fname);
  fprintf(stdout, "INFO: output file = %s\n", out_fname);

  // set the FFT size
  int samples_per_frame_per_pol = frame_dsize / (2 * (bits_per_sample / 4)); // number of complex samples per frame per pol
  fprintf(stdout, "INFO: samples per frame per pol is %d\n", samples_per_frame_per_pol);
  int fft_size = FFT_SIZE;
  int samples_per_channel = selected_samples / fft_size;

  fprintf(stdout, "INFO: channelising into %d channels, %d samples in each channel (number of FFTs)\n", fft_size, samples_per_channel);

  // allocate host memory buffers
  // int total_frames = total_samples / samples_per_frame_per_pol;
  char *input_frame = (char *)malloc(frame_dsize);
  fftw_complex *input_pol0 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * selected_samples);
  fftw_complex *input_pol1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * selected_samples);
  fftw_complex *output_pol0 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * selected_samples);
  fftw_complex *output_pol1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * selected_samples);
  double *input_double;
  fftw_plan plan_pol0 = fftw_plan_dft_1d(fft_size, input_pol0, output_pol0, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan plan_pol1 = fftw_plan_dft_1d(fft_size, input_pol1, output_pol1, FFTW_FORWARD, FFTW_MEASURE);
  float *power_both_pols = (float *)malloc(selected_samples * sizeof(float));

  // mask for signed vs unsigned data (used for 8 and 16 bit cases)
  int data_mask;
  if (format == 0) // signed data
    data_mask = 0;
  else                                      // unsigned data
    data_mask = 1 << (bits_per_sample - 1); // will be 0x80 for 8-bit data, 0x8000 for 16-bit data

  // read the input file and process all available frames
  size_t bytes_read;
  long int frames_read = 0;
  int selected_frame = 0;

  for (int frame = 0; frame < total_frames; frame++)
  {
    bytes_read = fread(input_frame, 1, frame_hsize, fp_in); // read the header bytes - which are ignored
    bytes_read = fread(input_frame, 1, frame_dsize, fp_in); // overwrite with the data

#if DEBUG
    if (frame == 0)
    {
      // print the first few bytes of the first frame
      fprintf(stdout, "First 16 bytes of first data frame:\n");
      for (int i = 0; i < 16; i++)
      {
        fprintf(stdout, "    byte[%d] = %d\n", i, (char)(input_frame[i] ^ data_mask));
      }
    }
#endif

    if (bytes_read != (size_t)frame_dsize)
    {
      fprintf(stdout, "INFO: cannot read %d bytes from input file, assuming reached EOF\n", frame_dsize);
    }
    else
    {
      frames_read++;
      // copy the input data to the two pol sample buffers IF we are a frame we want
      if (frame >= selected_start_frame && frame < selected_end_frame)
      {
        for (int i = 0; i < samples_per_frame_per_pol; i++)
        {
          input_double = (double *)(input_pol0 + (selected_frame * samples_per_frame_per_pol));
          input_double[2 * i] = (double)(char)(input_frame[4 * i + 0] ^ data_mask);
          input_double[2 * i + 1] = (double)(char)(input_frame[4 * i + 1] ^ data_mask);

          input_double = (double *)(input_pol1 + (selected_frame * samples_per_frame_per_pol));
          input_double[2 * i] = (double)(char)(input_frame[4 * i + 2] ^ data_mask);
          input_double[2 * i + 1] = (double)(char)(input_frame[4 * i + 3] ^ data_mask);
        }
        selected_frame++;
      }
    }
  }
  fprintf(stdout, "INFO: total number of frames read = %ld\n", frames_read);
  fprintf(stdout, "INFO: SELECTED frames read = %d\n", selected_frame);

  // free unneeded memory
  free(input_frame);

  // Close input file
  fclose(fp_in);

#if DEBUG
  // print the first samples of the the input file
  fprintf(stdout, "First complex voltage samples of the input file:\n");
  for (int i = 0; i < 10; i++)
  {
    input_double = (double *)input_pol0;
    fprintf(stdout, "    pol0 sample[%d] = %f   %f\n", i, input_double[2 * i], input_double[2 * i + 1]);
    input_double = (double *)input_pol1;
    fprintf(stdout, "    pol1 sample[%d] = %f   %f\n", i, input_double[2 * i], input_double[2 * i + 1]);
  }
#endif

  // FFT channelise each pol
  for (int ifft = 0; ifft < samples_per_channel; ifft++)
  {
    fftw_complex *in_ifft = input_pol0 + (ifft * fft_size);
    fftw_complex *out_ifft = output_pol0 + (ifft * fft_size);
    fftw_execute_dft(plan_pol0, in_ifft, out_ifft);
    in_ifft = input_pol1 + (ifft * fft_size);
    out_ifft = output_pol1 + (ifft * fft_size);
    fftw_execute_dft(plan_pol1, in_ifft, out_ifft);
  }
  fprintf(stdout, "INFO: %d FFTs completed for each pol\n", samples_per_channel);

  // destroy the FFTW stuff
  fftw_destroy_plan(plan_pol0);
  fftw_destroy_plan(plan_pol1);
  fftw_free(input_pol0);
  fftw_free(input_pol1);

#if DEBUG
  double *output_double;
  // print the first samples of the the channelised input file
  fprintf(stdout, "First complex voltage samples of the channelised input file:\n");
  for (int i = 0; i < 10; i++)
  {
    output_double = (double *)output_pol0;
    fprintf(stdout, "    pol0 channelised[%d] = %f   %f\n", i, output_double[2 * i], output_double[2 * i + 1]);
    output_double = (double *)output_pol1;
    fprintf(stdout, "    pol1 channelised[%d] = %f   %f\n", i, output_double[2 * i], output_double[2 * i + 1]);
  }
#endif

  // compute the power of both polarisations and sum them
  fprintf(stdout, "INFO: calculating the power of both polarisations and summing them\n");
  for (int i = 0; i < selected_samples; i++)
  {
    power_both_pols[i] = (float)(creal(output_pol0[i]) * creal(output_pol0[i]) + cimag(output_pol0[i]) * cimag(output_pol0[i]) +
                                 creal(output_pol1[i]) * creal(output_pol1[i]) + cimag(output_pol1[i]) * cimag(output_pol1[i]));
  }

  // Free final fftw variables
  fftw_free(output_pol0);
  fftw_free(output_pol1);

#if DEBUG
  // print the first powers
  fprintf(stdout, "First powers:\n");
  for (int i = 0; i < FFT_SIZE; i++) // one whole FFT size worth
  {
    fprintf(stdout, "    power both pols[%d] = %f\n", i, power_both_pols[i]);
  }
#endif

  // at this point we have Stokes I power samples in power_both_pols[] array, channelised

  fprintf(stdout, "INFO: Folding each fine channel\n");
  float *phase_bins = (float *)calloc(num_phase_bins, sizeof(float));         // bin values all initialised to zero
  int *num_values_per_bin = (int *)calloc(num_phase_bins, sizeof(int));       // bin counts all initialised to zero
  float *channel_dm_time_offsets = (float *)malloc(fft_size * sizeof(float)); // in seconds
  float fine_chan_width_MHz = (1.0E-6 * (float)SAMPLE_RATE) / (float)fft_size;
  float max_offset = -10.0; // set to a large negative number initially
  int channel;              // the channel number in ascending sky frequency order
  for (int c = 0; c < fft_size; c++)
  {
    if (c <= (fft_size / 2))
      channel = c;
    else
      channel = c - fft_size;
    float channel_freq_MHz = centre_frequency_MHz + ((float)channel * fine_chan_width_MHz);
    if (channel_freq_MHz != 0.0f)
      channel_dm_time_offsets[c] = 4.148808E3 * folding_dm / (channel_freq_MHz * channel_freq_MHz); // in seconds (4.148808E3 is dispersion constant MHz^2 cm^3 pc^-1 s)
    else
      channel_dm_time_offsets[c] = 0.0f;
    if (channel_dm_time_offsets[c] > max_offset)
      max_offset = channel_dm_time_offsets[c];
  }
  // adjust the offsets so that the largest is zero, and all others are negative
  for (int c = 0; c < fft_size; c++)
  {
    channel_dm_time_offsets[c] -= max_offset;
    if (channel_dm_time_offsets[c] > 0.0f) // just in case of rounding error
      channel_dm_time_offsets[c] = 0.0f;
  }
#if DEBUG // debug
  for (int c = 0; c < fft_size; c++)
  {
    fprintf(stdout, "channel %4d: DM time offset = %.6f s\n", c, channel_dm_time_offsets[c]);
  }
#endif

  double sample_period = (double)fft_size / SAMPLE_RATE; // time between successive samples in each channel

  // fold each channel separately, applying incoherent de-dispersion if requested
  int i, j;
  for (i = 0; i < fft_size; i++)
  {
    double start_time = phase_offset_time - channel_dm_time_offsets[i]; // in seconds, should not be negative (offsets are <= 0.0)

    // apply incoherent de-dispersion
    if (start_time < 0.0) // guarantee not negative (in case of rounding errors)
      start_time = 0.0;

    // fold each sample in this channel
    // - first sample is at start_time, increment by sample_period on each successive sample
    for (j = 0; j < samples_per_channel; j++)
    {
      double fractional_period = fmod(start_time, folding_period) / folding_period; // range is [0.0, 1.0)], with the value 1.0 possible but unlikely
      int phase_bin = (int)floor(fractional_period * (double)num_phase_bins);       // range is [0, num_phase_bins)], with the value num_phase_bins possible but unlikely
      if (phase_bin == num_phase_bins)
        phase_bin = num_phase_bins - 1; // ensures range is [0, num_phase_bins-1]
      num_values_per_bin[phase_bin]++;
      phase_bins[phase_bin] += power_both_pols[j * fft_size + i];

      // next sample
      start_time += sample_period;
    }
  }

  // free the  and DM time offsets array as we no longer need it
  free(channel_dm_time_offsets);
  free(power_both_pols);

  // scale the output profile by the number of values added to each bin
  for (i = 0; i < num_phase_bins; i++)
  {
    phase_bins[i] = phase_bins[i] / (float)num_values_per_bin[i];
  }

  if (normalise_output)
  {
    // find the maximum bin value
    float max_bin_value = -1.0f;
    for (i = 0; i < num_phase_bins; i++)
    {
      if (phase_bins[i] > max_bin_value)
        max_bin_value = phase_bins[i];
    }
    if (max_bin_value == 0.0)
      max_bin_value = 1.0; // don't normalise in zero case (avoids nan errors)

    // normalise all bins to the maximum value
    float normalisation_factor = 1.0f / max_bin_value;
    for (i = 0; i < num_phase_bins; i++)
    {
      phase_bins[i] *= normalisation_factor;
    }
  }

  // print all the bin values, and also save to the output file
  fprintf(stdout, "INFO: Pulse profile:\n");
  for (i = 0; i < num_phase_bins; i++)
  {
    if (log_output)
    {
      phase_bins[i] = 10.0f * log10f(phase_bins[i] + 1.0e-20f); // add a small value to avoid log10(0)
    }

    fprintf(stdout, "  bin[%3d]: %f   (num values = %d)\n", i, phase_bins[i], num_values_per_bin[i]);
    fprintf(fp_out, "%f\n", phase_bins[i]);
  }

  // free unneeded memory
  free(num_values_per_bin);

  // Do a plot
  double *x_axis = (double *)malloc(num_phase_bins * sizeof(double));
  for (i = 0; i < num_phase_bins; i++)
  {
    x_axis[i] = (double)i;
  }

  double *y_axis = (double *)malloc(num_phase_bins * sizeof(double));
  for (i = 0; i < num_phase_bins; i++)
  {
    y_axis[i] = (double)phase_bins[i];
  }

  // Free phase_bins after copying to y_axis
  free(phase_bins);

  char title[200];

  sprintf(title, "Folded Pulse Profile of %d ch: %03d beam: %02d (%d sec; %d to %d)", obsid, rec_channel, beam, selected_seconds, start_time, end_time);

  char plot_filename[200];
  sprintf(plot_filename, "pulse_profile_%d_ch%03d_beam%02d.png", obsid, rec_channel, beam);

  int ret = line_chart_png(plot_filename,
                           title,
                           "",
                           "Power",
                           x_axis,
                           y_axis,
                           num_phase_bins);

  // Free the x_axis and y_axis after plotting
  free(x_axis);
  free(y_axis);

  if (ret == 0)
    fprintf(stdout, "INFO: pulse profile plot written to %s\n", plot_filename);

  // Close out file
  fclose(fp_out);

  // final message
  fprintf(stdout, "INFO: all done, exiting\n\n");
  return EXIT_SUCCESS;
}
