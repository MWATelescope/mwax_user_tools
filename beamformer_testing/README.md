# fold_and_plot.sh

A bash script to automate pulsar folding and plotting from VDIF observation data using the `cirapulsarsandtransients/psr-analysis` Docker image.

## Overview

For a given set of beamformer observations that prodced VDIF and associated HDR files, the script will:

1. Run `dspsr` to fold each HDR file against a pulsar ephemeris (`.par` file), producing folded archive (`.ar`) files
2. Combine the resulting `.ar` files using `psradd` (skipped if only one HDR file is found)
3. Generate an integrated pulse profile plot using `pav`

## Requirements

- Docker
- Input directory containing `.vdif`, `.hdr`
- `.par` file to use for ephemeris. Get this by copying and pasting the "short" ephemeris from the [ATNF Pulsar Catalog](https://www.atnf.csiro.au/research/pulsar/psrcat/)
- The `cirapulsarsandtransients/psr-analysis` Docker image

## Usage

```bash
./fold_and_plot.sh \
  --data-dir <dir> \
  --output-dir <dir> \
  --start-obsid <id> \
  --end-obsid <id> \
  --par <local path to file> \
  --beam <num> \
  --chan <num> \
  [--threads <num>]
```

## Arguments

| Argument | Required | Description |
|---|---|---|
| `--data-dir` | Yes | Directory containing input VDIF, HDR, and PAR files |
| `--output-dir` | Yes | Directory to write output `.ar` and `.png` files |
| `--start-obsid` | Yes | 10-digit observation ID to start from (inclusive) |
| `--end-obsid` | Yes | 10-digit observation ID to end at (inclusive) |
| `--par` | Yes | PAR filename (filename, or full local path and filename) |
| `--beam` | Yes | Zero-padded 2 digit beam number (e.g. `01`, `02`, etc) |
| `--chan` | Yes | Zero-padded 3-digit receiver channel number (e.g. `091`, `123`, etc) |
| `--threads` | No | Number of threads for `dspsr` (default: `2`) |

## Output

All output files are written to `--output-dir`:

| File | Description |
|---|---|
| `<obsid>_ch<channel>_beam<beam>.ar` | Folded archive for each HDR file |
| `beam<beam>_<par>_combined.ar` | Combined archive across all observations |
| `beam<beam>_<par>_profile.png` | Integrated pulse profile plot |

## Example

```bash
./fold_and_plot.sh \
  --data-dir /data/observations \
  --output-dir /data/output \
  --start-obsid 1234567890 \
  --end-obsid 1234567899 \
  --par /home/myuser/J0835-4510.par \
  --beam 01 \
  --chan 091 \
  --threads 4
```

## Notes

- Only HDR files matching the specified channel, beam, and obsid range are processed
- The first HDR file is processed with an additional `-S 4` flag to `dspsr` to skip the first 4 seconds to account for QUACKTIME
- Output files are owned by the user running the script (not root)
- Run `./fold_and_plot.sh --help` for a quick usage summary

