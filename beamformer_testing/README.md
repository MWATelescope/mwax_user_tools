# run.sh

A wrapper around fold_and_plot.sh which allows running fold_and_plot across multiple obsids and channels.

## run.sh: Overview

This is just a wrapper that called fold_and_plot.sh 1-N times. It folds and then produces plots.

## run.sh: Usage

```bash
./run.sh --data-dir DIR --output-dir DIR --par PSR --start-obsid ID --end-obsid ID --beam BEAM --start-chan N --end-chan N [--bins N] [--threads N]
```

## run.sh: Arguments

| Argument | Required | Description |
|---|---|---|
| `--data-dir` | Yes | Path to the MWAX Beamformer VDIF files |
| `--output-dir` | Yes | Directory to write output `.ar` and `.png` files |
| `--start-obsid` | Yes | 10-digit observation ID to start from (inclusive) |
| `--end-obsid` | Yes | 10-digit observation ID to end at (inclusive) |
| `--par` | Yes | PAR filename (filename, or full local path and filename) |
| `--beam` | Yes | Zero-padded 2 digit beam number (e.g. `01`, `02`, etc) |
| `--start-chan` | Yes | Zero-padded 3-digit receiver channel number (e.g. `109`) |
| `--end-chan` | Yes | Zero-padded 3-digit receiver channel number (e.g. `132`) |
| `--threads` | No | Number of threads for `dspsr` (default: `4`) |
| `--bins` | No | Number of profile bins for dspsr (default: 64) |

# fold_and_plot.sh

A bash script to automate pulsar folding and plotting from VDIF observation data using the `cirapulsarsandtransients/psr-analysis` Docker image.
**NOTE** This script is run multiple times by the "run.sh" script, so normally you don't need to worry about it.

## fold_and_plot.sh: Overview

For a given set of beamformer observations that produced VDIF and associated HDR files, the script will:

1. Run `dspsr` to fold each HDR file against a pulsar ephemeris (`.par` file), producing folded archive (`.ar`) files
2. Combine the resulting `.ar` files using `psradd` (skipped if only one HDR file is found)
3. Generate an integrated pulse profile plot using `pav`

## fold_and_plot.sh: Requirements

- Docker
- Input directory containing `.vdif`, `.hdr`
- `.par` file to use for ephemeris. Get this by copying and pasting the "short" ephemeris from the [ATNF Pulsar Catalog](https://www.atnf.csiro.au/research/pulsar/psrcat/)
- The `cirapulsarsandtransients/psr-analysis` Docker image

## fold_and_plot.sh: Usage

```bash
./fold_and_plot.sh --data-dir <dir> --output-dir <dir> --start-obsid <id> --end-obsid <id> --par <file> --beam <num> --chan <num> [--threads <num>] [--bins <num>]
```

## fold_and_plot.sh: Arguments

| Argument | Required | Description |
|---|---|---|
| `--data-dir` | Yes | Directory containing input VDIF, HDR, and PAR files |
| `--output-dir` | Yes | Directory to write output `.ar` and `.png` files |
| `--start-obsid` | Yes | 10-digit observation ID to start from (inclusive) |
| `--end-obsid` | Yes | 10-digit observation ID to end at (inclusive) |
| `--par` | Yes | PAR filename (filename, or full local path and filename) |
| `--beam` | Yes | Zero-padded 2 digit beam number (e.g. `01`, `02`, etc) |
| `--chan` | Yes | Zero-padded 3-digit receiver channel number (e.g. `091`, `123`, etc) |
| `--threads` | No | Number of threads for `dspsr` (default: `4`) |
| `--bins` | No | Number of profile bins for dspsr (default: 64) |

## fold_and_plot.sh: Output

All output files are written to `--output-dir`:

| File | Description |
|---|---|
| `<obsid>_ch<channel>_beam<beam>.ar` | Folded archive for each HDR file |
| `<par>_<channel>_beam<beam>_combined.ar` | Combined archive across all observations |
| `<par>_<channel>_beam<beam>_combined_profile.png` | Integrated pulse profile plot |
| `<par>_<channel>_beam<beam>_combined_waterfall.png` | Waterfall plot |

## fold_and_plot.sh: Notes

- Only HDR files matching the specified channel, beam, and obsid range are processed
- The first HDR file is processed with an additional `-S 4` flag to `dspsr` to skip the first 4 seconds to account for QUACKTIME
- Output files are owned by the user running the script (not root)
- Run `./fold_and_plot.sh --help` for a quick usage summary

