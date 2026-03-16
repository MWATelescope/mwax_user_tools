# fold_vdif: VDIF Folding utility

## Building

```bash
gcc -o fold_vdif fold_vdif.c plot.c -lfftw3 -lm
```

## Example

This example will:
* Skip the first 4 seconds
* Specify the VDIF file has a 32 byte header in each frame
* Specify the VDIF file has a 8000 byte data block in each frame
* Specify the VDIF uses 8 bits per sample
* Speicfy the VDIF data is unsigned
* Speicfy 256 phase bins for folding
* Fold on a period of 1.2381 seconds
* Use 40.9 as the dispersion measure
* Move the period offset by 0.12 (shifts the pulse left or right)
* Specify central frequency of 168.96 MHz (based on the central frequency of the coarse channel)
* Output the pulse profile data as `pulse_profile.dat`
* Output the pulse profile plot as `pulse_profile_OBSID_chCCC_beamBB.png"`

```bash
./fold_vdif -s 4 -H 32 -d 8000 -b 8 -u -B 256 -P 1.2381 -D 40.9 -O 0.12 -f 168.96 -i /data/1455894016/1455894016_ch132_beam01.vdif
```

## Usage
```text
fold_vdif - fold a VDIF file at a specified period and DM

  Usage: fold_vdif <options>
         -s <int>          start time in seconds for folding [default 0 s]
         -e <int>          end time in seconds for folding [default -1 s]
         -H <int>          header size in bytes (0 for raw data, 32 for VDIF, 64 for CODIF) [default 32]
         -d <int>          data frame size in bytes [default 8000]
         -b <int>          bits per sample (real and imag) [default 8]
         -u                if this option is present, data is assumed to be unsigned (default is signed)
         -B <int>          number of phase bins for folding the data [default 128]
         -P <period>       period at which to fold the data (diagnostic) [default no folding]
         -D <num>          DM for incoherent de-dispersion [default no de-dispersion]
         -O <offset>       non-negative phase offset time in seconds (to place the pulse at a desired phase) [default 0.0 s]
         -f <freq>         centre frequency in MHz (for de-dispersion) [default 150.000000 MHz]
         -l                output the pulse profile in log scale (dB)
         -n                number of frames to process - minimum 1 [default all frames in the input file]
         -N                normalise the output pulse profile peak to 1.0 (0.0 for log output)
         -i                input filename
         -o                output filename [default pulse_profile.dat]
         -h                print this usage information
```