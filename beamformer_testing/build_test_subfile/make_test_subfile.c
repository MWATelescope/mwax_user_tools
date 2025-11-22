// A program to construct a test MWAX subfile from a set of signal files
// - taking header, signal and delay info from files specified on the command line

// The signal files are assumed to be raw complex samples of 8-bit signed integer
// with real and imaginary parts interleaved.  They are assumed to contain a full 8 seconds
// of data for an MWA sample rate of 1.28 Msample/s.

// The PSRDADA header, which is pre-pended to the output file, is taken from a text file
// specified on the command line.  If no header file is specified, a basic PSRDADA header
// is written to the output, just sufficient to allow the file to be loaded into a ring
// buffer using the dada_disk2db utility.
// (This header can later be overwitten using the update_MWA_header utility.)

// Build: gcc -o make_test_subfile make_test_subfile.c

// system include files
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>

#define PSRDADA_HEADER_SIZE 4096
#define FFT_SIZE 6400                // 5 ms, 10 FFTs per 50 ms block within a subfile
#define NUM_FFTS_PER_BLOCK 10
#define NUM_SUBBLOCKS 161            // 160 data sub-blocks for 8 second subfile, plus metadata in block 0
#define NUM_DELAYS 1600              // 1600 delays for 8 seconds of data, one per 5 ms (every FFT)
#define DELAY_ROW_BYTES 6456         // 56 bytes + 1600 * 4 bytes for the delays
#define NUM_INPUT_BITS 8             // twice this for a complex sample
#define DEFAULT_MWA_TILES 128
#define DEFAULT_SOURCES_FILENAME "subfile_source_files.txt"
#define DEFAULT_DELAYS_FILENAME "subfile_source_delays.txt"
#define DEFAULT_OUTPUT_FILENAME "test_subfile.sub"

void usage()
{
  fprintf(stdout, "\nmake_test_subfile - a program to construct a test MWAX subfile\n\n"
    "  Usage: make_test_subfile <options>\n"
    "         -T <int>          the number of dual-pol MWA tiles [default %d]\n"
    "         -H                header filename\n"
    "         -s <filename>     file containing filenames of input source signals for each tile/pol [default %s]\n"
    "         -D <filename>     file containing the delays that were applied generating each source [default %s]\n"
    "         -o                output subfile name [default %s]\n"
    "         -h                print this usage information\n\n",
    DEFAULT_MWA_TILES, DEFAULT_SOURCES_FILENAME, DEFAULT_DELAYS_FILENAME, DEFAULT_OUTPUT_FILENAME);
}


int main(int argc, char *argv[])
{
  FILE *fp_header;
  FILE *fp_sources;
  FILE *fp_delays;
  FILE *fp_out;

  // command line option parameters
  struct option options[] = {
    {"T",           required_argument, 0, 'T'},
    {"H",           required_argument, 0, 'H'},
    {"s",           required_argument, 0, 's'},
    {"D",           required_argument, 0, 'D'},
    {"o",           required_argument, 0, 'o'},
    {"h",           no_argument,       0, 'h'},
    {0, 0, 0, 0}
  };

  int opt, ss;
  char header_fname[100];
  header_fname[0]= '\0';  // empty string
  char sources_fname[100] = DEFAULT_SOURCES_FILENAME;
  char delays_fname[100] = DEFAULT_DELAYS_FILENAME;
  char out_fname[100] = DEFAULT_OUTPUT_FILENAME;

struct block_0_path_metadata_sans_delays
{
  uint16_t rf_input;
  int16_t ws_delay_applied;    // the whole-sample delay for this signal path, for the entire sub-observation - will often be negative
  double start_total_delay;    // the start, middle and delays, from which intermediate delays can be calculated by polynomial interpolation
  double middle_total_delay;
  double end_total_delay;
  double initial_delay;        // initial residual delay at centre of first 5 ms timestep
  double delta_delay;          // increment between timesteps
  double delta_delta_delay;    // increment in the increment between timesteps
  int16_t num_pointings;       // may be multiple pointings with beamforming, the first pointing is the correlation pointing centre
  int16_t reserved;
 };

  int tiles = DEFAULT_MWA_TILES;
  int fft_size = FFT_SIZE;
  int num_input_bits = NUM_INPUT_BITS;

  int bits_per_sample = 2 * num_input_bits;  // complex samples

  // parse command line options
  while (1)
  {
    opt = getopt_long_only(argc, argv, "T:H:s:D:o:h", options, NULL);
    
    if (opt == EOF) break;
    
    switch (opt)
    {
      case 'T':
        tiles = atoi(optarg);
        break;

      case 'H':
        ss = sscanf(optarg, "%s", header_fname);
        if (ss!=1)
        {
          fprintf(stderr, "ERROR: bad -file option %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 's':
        ss = sscanf(optarg, "%s", sources_fname);
        if (ss!=1)
        {
          fprintf(stderr, "ERROR: bad -file option %s\n", optarg);
          exit(EXIT_FAILURE);
        }
        break;

      case 'D':
        ss = sscanf(optarg, "%s", delays_fname);
        if (ss!=1)
        {
          fprintf(stderr, "ERROR: bad -file option %s\n", optarg);
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

  // check for a header file, and if there is one specified, try to open it
  if (header_fname[0] == '\0')
  {
    fprintf(stderr, "WARNING: no header file specified (with the -H option), basic default header will be written\n");
  }
  else
  {
    fp_header = fopen(header_fname, "r");
    if (fp_header == NULL)
    {
      fprintf(stderr, "ERROR: cannot open header file %s\n", header_fname);
      exit(EXIT_FAILURE);
    }
  }
  // open the sources filename file
  fp_sources = fopen(sources_fname, "r");
  if (fp_sources == NULL)
  {
    fprintf(stderr, "ERROR: cannot open source signal filenames file %s\n", sources_fname);
      exit(EXIT_FAILURE);
  }
  // open the delays filename file
  fp_delays = fopen(delays_fname, "r");
  if (fp_delays == NULL)
  {
    fprintf(stderr, "ERROR: cannot open source delays file %s\n", delays_fname);
    exit(EXIT_FAILURE);
  }
  // open the output file
  fp_out = fopen(out_fname, "w");
  if (fp_out == NULL)
  {
    fprintf(stderr, "ERROR: cannot open output file %s\n", out_fname);
    exit(EXIT_FAILURE);
  }

  int num_sources = tiles * 2;  // 2 pols per tile

  // try to open all the specified signal files
  // allocate an array of file pointers of size tiles
  FILE **fp_sources_array = malloc(num_sources * sizeof(FILE *));
  char **source_fname = malloc(num_sources * sizeof(char *));
  for (int i = 0; i < num_sources; i++)
  {
    fp_sources_array[i] = NULL;  // initialise to NULL
    source_fname[i] = malloc(100 * sizeof(char));  // allocate space for each filename
  }
  for (int i = 0; i < tiles; i++)
  {
    if (fscanf(fp_sources, "%s", source_fname[2*i]) != 1)
    {
      fprintf(stderr, "ERROR: cannot read pol0 signal filename for tile %d\n", i);
      exit(EXIT_FAILURE);
    }
    fp_sources_array[2*i] = fopen(source_fname[2*i], "r");
    if (fp_sources_array[2*i] == NULL)
    {
      fprintf(stderr, "ERROR: cannot open pol0 signal file %s for tile %d\n", source_fname[2*i], i);
      exit(EXIT_FAILURE);
    }
    if (fscanf(fp_sources, "%s", source_fname[2*i + 1]) != 1)
    {
      fprintf(stderr, "ERROR: cannot read pol1 signal filename for tile %d\n", i);
      exit(EXIT_FAILURE);
    }
    fp_sources_array[2*i + 1] = fopen(source_fname[2*i + 1], "r");
    if (fp_sources_array[2*i + 1] == NULL)
    {
      fprintf(stderr, "ERROR: cannot open pol1 signal file %s for tile %d\n", source_fname[2*i + 1], i);
      exit(EXIT_FAILURE);
    }
  }

  if (header_fname[0] != '\0')
    fprintf(stdout, "INFO: input header file = %s\n", header_fname);
  for (int i = 0; i < tiles; i++)
  {
    fprintf(stdout, "INFO: input signal source file for tile %d, pol0 = %s\n", i, source_fname[2*i]);
    fprintf(stdout, "INFO: input signal source file for tile %d, pol1 = %s\n", i, source_fname[2*i + 1]);
  }
  fprintf(stdout, "INFO: input delays file = %s\n", delays_fname);
  fprintf(stdout, "INFO: output file = %s\n", out_fname);

  // set up the header in the output file
  char header[PSRDADA_HEADER_SIZE];
  int char_count = 0;
  if (header_fname[0] == '\0')
  {
    // write a basic PSRDADA header to the header array
    snprintf(header, PSRDADA_HEADER_SIZE, "HDR_SIZE %d\nNBIT %d\nUTC_START 2025-07-01-00:00:00\nOBS_OFFSET 0\n", PSRDADA_HEADER_SIZE, bits_per_sample);
    char_count = strlen(header);
  }
  else
  {
    // copy the header file to the header array, copying one char at a time until EOF or PSRDADA_HEADER_SIZE is reached
    while (char_count < PSRDADA_HEADER_SIZE && (header[char_count] = fgetc(fp_header)) != EOF) char_count++;
  }
  if (char_count < PSRDADA_HEADER_SIZE)
  {
    // fill the rest of the header with spaces
    memset(header + char_count, ' ', PSRDADA_HEADER_SIZE - char_count);
  }

  // write the header to the output file
  if (fwrite(header, 1, PSRDADA_HEADER_SIZE, fp_out) != PSRDADA_HEADER_SIZE)
  {
    fprintf(stderr, "ERROR: failed to write header to output file %s\n", out_fname);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout, "INFO: header written to output file %s\n", out_fname);

  // allocate memory for a buffer to hold all 161 sub-blocks of the subfile - initialised to zero
  int num_subblocks = NUM_SUBBLOCKS;
  // 64000 samples per sub-block, 50 ms at 1.28 Msample/s, 10 FFTs of 6400 samples each
  int subblock_row_size = NUM_FFTS_PER_BLOCK * fft_size * (bits_per_sample / 8);  // size in bytes
  long int subblock_total_size = subblock_row_size * num_sources;
  size_t total_size = num_subblocks * subblock_total_size;

  fprintf(stdout, "INFO: size of data portion of subfile = %zu bytes, size including header = %zu\n", total_size, total_size + PSRDADA_HEADER_SIZE);

  int8_t *subfile_buffer = calloc(total_size, sizeof(int8_t));
  if (subfile_buffer == NULL)
  {
    fprintf(stderr, "ERROR: failed to allocate memory for subfile buffer\n");
    exit(EXIT_FAILURE);
  }

  // loop through the input files, populating the data portion of the subfile buffer (note sub-block index starting at 1)
  for (int i = 1; i < num_subblocks; i++)
  {
    // for each sub-block, read the data from each source file
    for (int j = 0; j < num_sources; j++)
    {
      size_t bytes_read = fread(subfile_buffer + (i * subblock_total_size) + (j * subblock_row_size), 1, subblock_row_size, fp_sources_array[j]);
      if (bytes_read != subblock_row_size)
      {
        fprintf(stderr, "ERROR: failed to read enough data from source file %s for sub-block %d source %d: expected %d bytes, got %zu\n",
                source_fname[j], i, j, subblock_row_size, bytes_read);
        exit(EXIT_FAILURE);
      }
    }
  }

  fprintf(stdout, "INFO: data sub-blocks populated in the subfile buffer\n");

  // set up the metadata for the first sub-block (block 0)
  // set to all zeros for now
  memset(subfile_buffer, 0, subblock_total_size);

  // write the delay values from the delays file into block0
  #define SME_OFFSET 4           // where the start/middle/end delays start in each delay row
  #define DELAY_TABLE_OFFSET 56  // where the delay table values start in each delay row
  int delay_row_float_stride = DELAY_ROW_BYTES / sizeof(float);    // number of floats in a delay row
  int delay_row_double_stride = DELAY_ROW_BYTES / sizeof(double);  // number of doubles in a delay row
  double *sme_start = (double *)(subfile_buffer + SME_OFFSET);
  float *delays_start = (float *)(subfile_buffer + DELAY_TABLE_OFFSET);
  for (int i = 0; i < num_sources; i++)
  {
    float delay;
    if (fscanf(fp_delays, "%f", &delay) != 1)
    {
      fprintf(stderr, "ERROR: failed to read delay value for source %d\n", i);
      exit(EXIT_FAILURE);
    }
    // print the delay value
    fprintf(stdout, "INFO: delay for source %d = %f samples\n", i, delay);
    // set start/middle/end delays to the same value
    sme_start[i * delay_row_double_stride] = (double)delay;
    sme_start[i * delay_row_double_stride + 1] = (double)delay;
    sme_start[i * delay_row_double_stride + 2] = (double)delay;
    // populate the delay table with this value for all 1600 entries
    for (int k = 0; k < NUM_DELAYS; k++)
    {
      delays_start[i * delay_row_float_stride + k] = delay;
    }
  }

  // copy all sub-blocks of the buffer to the output file
  if (fwrite(subfile_buffer, 1, total_size, fp_out) != total_size)
  {
    fprintf(stderr, "ERROR: failed to write subfile buffer to output file %s\n", out_fname);
    exit(EXIT_FAILURE);
  }

  fprintf(stdout, "INFO: metadata and data sub-blocks written to output file %s\n", out_fname);

  if (header_fname[0] != '\0') fclose(fp_header);
  fclose(fp_sources);
  fclose(fp_delays);
  fclose(fp_out);
  for (int i = 0; i < num_sources; i++)
  {
    fclose(fp_sources_array[i]);
  }

  // free allocated memory
  for (int i = 0; i < num_sources; i++)
  {
    free(source_fname[i]);
  }
  free(source_fname);
  free(fp_sources_array);
  free(subfile_buffer);

  // final message
  fprintf(stdout, "INFO: all done, exiting\n\n");
  return EXIT_SUCCESS;
}
