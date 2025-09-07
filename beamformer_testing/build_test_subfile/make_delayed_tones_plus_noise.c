// A program to generate a test signal consisting of:
// - up to 3 tones, generated in the frequency domain, with specified FFT bins,
// - a delay applied to the signal, specified as a fraction of one sample,
// - added Guassian noise, with a specified noise amplitude.
//
// Each output frame has a header of selectable length.  Choose length 0 for contiguous raw data frames,
// 32 for VDIF frames, 64 for CODIF frames, or whatever length you want.
// NOTE: the headers are unpopulated (all zeros).
//
// The frame size and FFT size are chosen together: the FFT size is always the number of complex samples per frame.
//
// If a second tone is specified, it is always included at half the amplitude of the first tone (but this can be changed by hand).
// If a third tone is specified, it is included at one quarter the amplitude of the first tone (but this can be changed by hand).
//
// A command line parameter sets the delay of the signal wrt non-delayed, specified as a fraction of one sample.
//
// Build with: gcc -o make_delayed_tones_plus_noise make_delayed_tones_plus_noise.c -lfftw3 -lm
//
// Ian Morrison July 2025

// system include files
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <fftw3.h>

#define DEFAULT_HEADER_SIZE 0             // default is contiguous raw data frames
#define DEFAULT_DATA_FRAME_SIZE 8192      // 8192 bytes
#define DEFAULT_NUM_INPUT_BITS 8          // twice this for a complex sample
#define DEFAULT_TONE_BIN 0                // default FFT bin for the first tone
#define SAMPLE_RATE_MHZ 1280000.0         // 1.28 MHz
#define DEFAULT_OUTPUT_FILENAME "output_samples.cdf"

void usage()
{
  fprintf(stdout, "\nmake_delayed_tones_plus_noise - a tool to generate a noisy multi-tone test signal in VDIF, CODIF or raw data frame format\n\n"
    "  Usage: nmake_delayed_tones_plus_noise <options>\n"
    "         -H <int>          header size in bytes (e.g. 0 for raw data, 32 for VDIF, 64 for CODIF) [default %d]\n"
    "         -d <int>          data frame size in bytes [default %d]\n"
    "         -b <int>          bits per sample (real and imag) [default %d]\n"
    "         -1                FFT bin for first tone [default %d]\n"
    "         -2                FFT bin for second tone (if not specified, only one tone is generated)\n"
    "         -3                FFT bin for third tone (if not specified, only one or two tones are generated)\n"
    "         -D <float>        delay of the signal in samples (fractional) [default 0.0]\n"
    "         -l                length of the output file in seconds [default 8 seconds]\n"
    "         -n <float>        amplitude of the noise to add to the signal [default 0.0]\n"
    "         -s <int>          seed (positive integer) for the random number generator [default 1]\n"
    "         -o                output filename [default %s]\n"
    "         -h                print this usage information\n\n",
    DEFAULT_HEADER_SIZE, DEFAULT_DATA_FRAME_SIZE, DEFAULT_NUM_INPUT_BITS, DEFAULT_TONE_BIN, DEFAULT_OUTPUT_FILENAME);
}



int main(int argc, char *argv[])
{
  FILE *fp_out;

  // command line option parameters
  struct option options[] = {
    {"H",           required_argument, 0, 'H'},
    {"d",           required_argument, 0, 'd'},
    {"b",           required_argument, 0, 'b'},
    {"1",           required_argument, 0, '1'},
    {"2",           required_argument, 0, '2'},
    {"3",           required_argument, 0, '3'},
    {"D",           required_argument, 0, 'D'},
    {"l",           required_argument, 0, 'l'},
    {"n",           required_argument, 0, 'n'},
    {"s",           required_argument, 0, 's'},
    {"o",           required_argument, 0, 'o'},
    {"h",           no_argument,       0, 'h'},
    {0, 0, 0, 0}
  };

  int opt, ss;
  char out_fname[100] = DEFAULT_OUTPUT_FILENAME;
  double noise_amplitude = 0.0; // default noise amplitude (noise-less)
  int seed = 1;                 // default seed for random number generator
  int tone_bin[3] = { DEFAULT_TONE_BIN, -1, -1 };  // first tone bin is set to default, others unset
  double tone_amplitude[3] = { 1.0, 0.5, 0.25 };
  int frame_hsize = DEFAULT_HEADER_SIZE;
  int frame_dsize = DEFAULT_DATA_FRAME_SIZE;
  int bits_per_sample = DEFAULT_NUM_INPUT_BITS;  // 2x this for full complex sample
  double delay_fraction = 0.0;  // default delay is zero
  float length_seconds = 8.0;   // default length of the output file in seconds

  // parse command line options
  while (1)
  {
    opt = getopt_long_only(argc, argv, "H:d:b:1:2:3:D:l:n:s:o:h", options, NULL);
    
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
        break;

      case '1':
        tone_bin[0] = atoi(optarg);
        break;

      case '2':
        tone_bin[1] = atoi(optarg);
        break;

      case '3':
        tone_bin[2] = atoi(optarg);
        break;

      case 'D':
        delay_fraction = atof(optarg);
        if ((delay_fraction < 0.0) || (delay_fraction > 1.0))
        {
          fprintf(stderr, "ERROR: delay fraction must be in the range [0.0, 1.0), got %f\n", delay_fraction);
          exit(EXIT_FAILURE);
        }
        break;

      case 'l':
        length_seconds = atof(optarg);
        if (length_seconds <= 0.0)
        {
          fprintf(stderr, "ERROR: length must be a positive number, got %f\n", length_seconds);
          exit(EXIT_FAILURE);
        }
        break;

      case 'n':
        noise_amplitude = atof(optarg);
        if (noise_amplitude < 0.0)
        {
          fprintf(stderr, "ERROR: noise amplitude must be non-negative, got %f\n", noise_amplitude);
          exit(EXIT_FAILURE);
        }
        break;

      case 's':
        seed = atoi(optarg);
        if (seed < 1)
        {
          fprintf(stderr, "ERROR: seed must be a positive integer, got %d\n", seed);
          exit(EXIT_FAILURE);
        }
        break;

      case 'o':
        ss = sscanf(optarg, "%s", out_fname);
        if (ss!=1)
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

  int bits_per_frame = 8 * frame_dsize;        // frame_dsize is in bytes
  if (bits_per_frame % (2*bits_per_sample))
  {
    fprintf(stderr, "ERROR: cannot fit an integer number of complex samples (each 2 x %d bits) into a single frame of size %d bits, exiting\n", bits_per_sample, bits_per_frame);
    exit(EXIT_FAILURE);
  }

  // open the output file
  fp_out = fopen(out_fname, "w");
  if (fp_out == NULL)
  {
    fprintf(stderr, "ERROR: cannot open output file %s\n", out_fname);
    exit(EXIT_FAILURE);
  }
  fprintf(stdout, "INFO: output file = %s\n", out_fname);

  // set the FFT size
  int samples_per_frame = frame_dsize / (bits_per_sample/4);  // number of complex samples per frame
  int fft_size = samples_per_frame;
  fprintf(stdout, "INFO: number of samples per frame = FFT size = %d\n", fft_size);

  // allocate host memory buffers
  fftw_complex *data_frame_complex_double = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*fft_size);
  void *out_frame = malloc(frame_dsize);
  char *out_frame_8bit = (char *)out_frame;
  short *out_frame_16bit = (short *)out_frame;

  // make an empty header
  char *header = (char *)calloc(frame_hsize, sizeof(char));

  // create the inverse FFT plan
  fftw_plan plan = fftw_plan_dft_1d(fft_size, data_frame_complex_double, data_frame_complex_double, FFTW_BACKWARD, FFTW_MEASURE);

  // initialise the data frame with complex phasors of the specified amplitude and random phase
  // (bins with tones will be overwritten)
  srand(seed);  // seed the random number generator
  fprintf(stdout, "INFO: generating noise with amplitude %f and seed %d\n", noise_amplitude, seed);
  for (int i=0; i<fft_size; i++)
  {
    double random_phase = ((double)rand() / (double)RAND_MAX) * 2.0 * M_PI; // random phase in [0, 2pi]
    data_frame_complex_double[i][0] = noise_amplitude * cos(random_phase);
    data_frame_complex_double[i][1] = noise_amplitude * sin(random_phase);
  }

  // count the number of tones set, checking they're in range
  int num_tones = 0;
  for (int i=0; i<3; i++)
  {
    if (tone_bin[i] == -1) continue;  // skip unset tone bins

    if ((tone_bin[i] >= 0) && (tone_bin[i] < fft_size))
    {
      fprintf(stdout, "INFO: tone bin %d = %d\n", (i+1), tone_bin[i]);
      num_tones++;
    }
    else
    {
      fprintf(stderr, "ERROR: tone bin %d (%d) is out of range [0, %d]\n", (i+1), tone_bin[i], (fft_size - 1));
      goto cleanup;
    }
  }

  // insert the specified tones at their specified amplitudes and phases
  // (negative phase shifts correspond to positive delays, i.e. shifting to later in time)
  double phase_shift_per_chan = -2.0 * M_PI * delay_fraction / (double)fft_size;
  fprintf(stdout, "INFO: applying a delay of %f samples to the signal\n", delay_fraction);
  for (int i=0; i<num_tones; i++)
  {
    double phase_shift_for_this_bin;
    if (tone_bin[i] <= (fft_size/2))
      phase_shift_for_this_bin = phase_shift_per_chan * (double)tone_bin[i];
    else
      phase_shift_for_this_bin = phase_shift_per_chan * (double)(tone_bin[i] - fft_size);
    data_frame_complex_double[tone_bin[i]][0] = tone_amplitude[i] * cos(phase_shift_for_this_bin);
    data_frame_complex_double[tone_bin[i]][1] = tone_amplitude[i] * sin(phase_shift_for_this_bin);
    fprintf(stdout, "INFO: tone bin %d phase shifted by %f radians = %f degrees\n", tone_bin[i], phase_shift_for_this_bin, (phase_shift_for_this_bin * (180.0/M_PI)));
  }

  // execute the FFT plan, which generates the time-domain signal
  fftw_execute(plan);

  // find the maximum real or imag value in the time-domain signal
  double max_value = 0.0;
  for (int i=0; i<fft_size; i++)
  {
    if (fabs(data_frame_complex_double[i][0]) > max_value)
      max_value = fabs(data_frame_complex_double[i][0]);
    if (fabs(data_frame_complex_double[i][1]) > max_value)
      max_value = fabs(data_frame_complex_double[i][1]);
  }

  // calculate the number of frames to generate
  long int num_frames = (long int)((length_seconds * SAMPLE_RATE_MHZ) / (float)samples_per_frame);

  // convert the data frame to nbits per sample format
  if (bits_per_sample == 8)
  {
    // convert to 8-bit samples
    for (int i=0; i<samples_per_frame; i++)
    {
      out_frame_8bit[2*i] = (char)((data_frame_complex_double[i][0]/max_value) * 127.0);   // real part
      out_frame_8bit[2*i+1] = (char)((data_frame_complex_double[i][1]/max_value) * 127.0); // imaginary part
    }
  }
  else if (bits_per_sample == 16)
  {
    // convert to 16-bit samples
    for (int i=0; i<samples_per_frame; i++)
    {
      out_frame_16bit[2*i] = (short)((data_frame_complex_double[i][0]/max_value) * 32767.0);   // real part
      out_frame_16bit[2*i+1] = (short)((data_frame_complex_double[i][1]/max_value) * 32767.0); // imaginary part
    }
  }
  else
  {
    fprintf(stderr, "ERROR: unsupported bits per sample %d, only 8 or 16 are supported\n", bits_per_sample);
    goto cleanup;
  }
  
  // write the header and data to the output file
  for (int i=0; i<num_frames; i++)
  {
    fwrite(header, frame_hsize, 1, fp_out);
    fwrite(out_frame, frame_dsize, 1, fp_out);
  }
  fprintf(stdout, "INFO: generated %ld frames of data in the output file\n", num_frames);

#if 1
  // print the first samples of the last frame of the output
  fprintf(stdout, "First complex voltage samples of the generated signal:\n");
  for (int i=0; i<10; i++)
  {
    if (bits_per_sample == 8)
      fprintf(stdout, "    sample[%d] = %d   %d\n", i, out_frame_8bit[2*i], out_frame_8bit[2*i+1]);
    else if (bits_per_sample == 16)
      fprintf(stdout, "    sample[%d] = %d   %d\n", i, out_frame_16bit[2*i], out_frame_16bit[2*i+1]);
    else
    {
      fprintf(stderr, "ERROR: unsupported bits per sample %d, must be 8 or 16\n", bits_per_sample);
      break;
    }
  }
#endif

  cleanup:

  fclose(fp_out);
  free(header);
  free(out_frame);

  // destroy the FFTW stuff
  fftw_destroy_plan(plan);
  fftw_free(data_frame_complex_double);

  // final message
  fprintf(stdout, "INFO: all done, exiting\n\n");
  return EXIT_SUCCESS;
}
