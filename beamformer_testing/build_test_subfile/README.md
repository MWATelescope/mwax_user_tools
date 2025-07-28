# Build Test Subfile

## build_test_subfile.sh

This is a bash script for constructing subfiles for beamformer testing.

To run:
```
$ ./build_test_subfile.sh
```

Edit the parameters at the top of the script to set the number of tiles and signal specifications.

It calls the C program "make_delayed_tones_plus_noise" multiple times, creating multiple versions of the same
signal, each with a different fractional sample delay.

It then calls another C program "make_test_subfile" to construct the subfile from the constituent signal files,
and populate the delays array in block 0 of the subfile.  A PSRDADA header is copied from a file specified
on the command line.

## make_delayed_tones_plus_noise.c

```
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

Usage: nmake_delayed_tones_plus_noise <options>
         -H <int>          header size in bytes (e.g. 0 for raw data, 32 for VDIF, 64 for CODIF) [default 0]
         -d <int>          data frame size in bytes [default 8192]
         -b <int>          bits per sample (real and imag) [default 8]
         -1                FFT bin for first tone [default 0]
         -2                FFT bin for second tone (if not specified, only one tone is generated)
         -3                FFT bin for third tone (if not specified, only one or two tones are generated)
         -D <float>        delay of the signal in samples (fractional) [default 0.0]
         -l                length of the output file in seconds [default 8 seconds]
         -n <float>        amplitude of the noise to add to the signal [default 0.0]
         -s <int>          seed (positive integer) for the random number generator [default 1]
         -o                output filename [default output_samples.cdf]
         -h                print this usage information
```

## make_test_subfile.c

```
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

Usage: make_test_subfile <options>
         -T <int>          the number of dual-pol MWA tiles [default 128]
         -H                header filename
         -s <filename>     file containing filenames of input source signals for each tile/pol [default subfile_source_files.txt]
         -D <filename>     file containing the delays that were applied generating each source [default subfile_source_delays.txt]
         -o                output subfile name [default test_subfile.sub]
         -h                print this usage information
```


## Compile both C programs:

```
$ make
```
