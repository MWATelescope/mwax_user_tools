use anyhow::{Context, Result, bail, ensure};
use fitsio::FitsFile;
use ndarray::ArrayD;
use num::Complex;
use std::env::args;

const ALLOWED_DIFFERENT: [&'static str; 5] = ["CORR_VER", "U2S_VER", "CBF_VER", "DB2F_VER", "TIME"];
const EXPECTED_IDENTICAL: [&'static str; 15] = [
    "SIMPLE", "BITPIX", "NAXIS", "MARKER", "MILLITIM", "PROJID", "OBSID", "FINECHAN", "NFINECHS",
    "INTTIME", "NINPUTS", "CORRHOST", "CORRCHAN", "MC_IP", "MC_PORT",
];

// fn wat(n1:&str) -> Result<()> {
//     let mut f1 = FitsFile::open(n1)?;
//     let mut f2 = FitsFile::open(n1)?;
//
//     for h1 in f2.iter() {
//         if let Ok(ninputs) = h1.read_key::<i64>(&mut f1, "NINPUTS") {
//             println!("number of inputs = {ninputs}");
//         }
//     }
//     Ok(())
// }

fn tri(x: usize) -> usize {
    x * (x + 1) / 2
}

fn cmp(n1: &str, n2: &str) -> Result<()> {
    let mut f1 = FitsFile::open(n1)?;
    let mut f2 = FitsFile::open(n2)?;
    let mut _failed = false;
    let mut identical = true;

    println!("comparing {}", n1);
    println!("      and {}", n2);

    let num_tiles;
    let int_time;
    let fine_channels;
    {
        let h1 = f1.hdu(0)?;
        let h2 = f2.hdu(0)?;
        for k in ALLOWED_DIFFERENT.iter() {
            let v1 = h1
                .read_key::<String>(&mut f1, k)
                .with_context(|| format!("no {k} in {n1} header"))?;
            let v2 = h2
                .read_key::<String>(&mut f2, k)
                .with_context(|| format!("no {k} in {n2} header"))?;
            println!("{:>8}: {:>8} <=> {:8}", k, v1, v2);
        }
        for k in EXPECTED_IDENTICAL.iter() {
            let v1 = h1
                .read_key::<String>(&mut f1, k)
                .with_context(|| format!("no {k} in {n1} header"))?;
            let v2 = h2
                .read_key::<String>(&mut f2, k)
                .with_context(|| format!("no {k} in {n2} header"))?;
            ensure!(v1 == v2, "{}: {:} != {:}", k, v1, v2);
            println!("{k}: {v1}")
        }
        num_tiles = h1.read_key::<i64>(&mut f1, "NINPUTS")? as usize / 2;
        int_time = h1.read_key::<f64>(&mut f1, "INTTIME")?;
        fine_channels = h1.read_key::<i64>(&mut f1, "NFINECHS")? as usize;
    }
    println!("number of tiles = {num_tiles}");
    println!("integration time = {int_time}");
    println!("fine_channels = {fine_channels}");
    let int_count = (8.0f64 / int_time) as usize;
    println!("integration count = {int_count}");

    let mut chan_errs: Vec<f32> = Vec::with_capacity(fine_channels);

    const EXTENSION_KEYWORDS: [&'static str; 3] = ["MARKER", "TIME", "MILLITIM"];

    let mut hdu_index = 1;
    loop {
        chan_errs.clear();
        chan_errs.resize(fine_channels, 0.0);
        let h1 = match f1.hdu(hdu_index) {
            Ok(h) => h,
            Err(_) => {
                if f2.hdu(hdu_index).is_ok() {
                    bail!(&"too many HDUs in file 2");
                }
                break; // Reached end of file 1 naturally, break loop safely
            }
        };
        let h2 = f2.hdu(hdu_index).context(format!(
            "File 2 is missing extension HDU index {hdu_index}!"
        ))?;

        for &k in EXTENSION_KEYWORDS.iter() {
            match (
                h1.read_key::<String>(&mut f1, k),
                h2.read_key::<String>(&mut f2, k),
            ) {
                (Err(e), Err(_)) => {
                    println!("both files, hdu {hdu_index} is missing {k} {e:?}")
                }
                (Ok(_), Err(e)) => {
                    println!("{n2} hdu {hdu_index} is missing {k} {e:?}")
                }
                (Err(e), Ok(_)) => {
                    println!("{n1} hdu {hdu_index} is missing {k} {e:?}")
                }
                (Ok(v1), Ok(v2)) => {
                    ensure!(
                        v1 == v2,
                        "Metadata Mismatch in HDU {hdu_index}! Key '{k}': reference='{v1}' != target='{v2}'"
                    );
                    print!(" {k}=>{v1:3}")
                }
            }
        }

        let d1: ArrayD<f32> = h1.read_image(&mut f1).unwrap();
        let d2: ArrayD<f32> = h2.read_image(&mut f2).unwrap();

        let nvis = d1.shape()[0];
        let vislen = d1.shape()[1];

        // 1. DYNAMIC EDGE-CASE ROUTING: Validate dimensions based on extension layer type
        if hdu_index % 2 == 1 {
            // ODD HDU INDEX == Visibilities Image Matrix
            ensure!(
                vislen == fine_channels * 8,
                "wrong visibility length in hdu {hdu_index} (expected {} got {vislen})",
                fine_channels * 8
            );
            ensure!(
                nvis == tri(num_tiles),
                "wrong visibility count in hdu {hdu_index} (expected {} got {nvis})",
                tri(num_tiles)
            );
        } else {
            // EVEN HDU INDEX == Weights Matrix
            ensure!(
                vislen == 4,
                "wrong weights length in hdu {hdu_index} (expected 4 got {vislen})"
            );
        }

        let mut m = 0.0;
        let mut max = 0.0f32;
        for visep in 0..nvis {
            for j in 0..vislen / 2 {
                let c1 = Complex::new(d1[[visep, j * 2 + 0]], d1[[visep, j * 2 + 1]]);
                let c2 = Complex::new(d2[[visep, j * 2 + 0]], d2[[visep, j * 2 + 1]]);
                let err = 2.0 * (c1 - c2).norm() / (c1.norm() + c2.norm() + 1e-38);
                m += err * err;
                max = max.max(err);
                chan_errs[j / 4] += err * err / ((nvis * vislen / 2) as f32);
            }
            let rms = (m / vislen as f32).sqrt();
            if rms < 0.0 {
                println!(
                    " hdu {hdu_index}: rms relative error = {rms:9.7}, max relative error = {max:9.7}"
                );
            }
            if max > 0.0001f32 || rms > 0.0001f32 {
                _failed = true;
            }
        }
        if max > 0.0f32 {
            identical = false;
        }
        let ner = chan_errs.len();
        let period = (hdu_index-1) / 2;
        let d = (hdu_index-1) % 2;
        println!(
            "  Period {period} datablock {d} fine channel mean squared errors: {:?}",
            &chan_errs[..ner.min(20)]
        );

        hdu_index += 1; // Advance to next timestep extension block smoothly
    }
    println!(
        "Successfully cross-checked {} total integration extensions.",
        hdu_index - 1
    );

    // ensure!(!_failed, "differences are too high\n");
    if identical {
        println!("inputs are identical!");
    }
    println!();

    Ok(())
}
fn usage() {
    println!("Usage:\n\tfitscmp <f1.fits> <f2.fits>");
}
fn main() -> Result<()> {
    unsafe {
        std::env::set_var("RUST_LIB_BACKTRACE", "0"); // bodge to kill the backtrace on mismatches
    }
    let args = args().collect::<Vec<String>>();
    if args.len() != 3 {
        usage();
        bail!("incorrect number of arguments");
    }
    //let f1 = "/Users/shrydar/cira/fits/ori/1436174864_20250710092726_ch070_000.fits";
    //let f2 = "/Users/shrydar/cira/fits/new/1436174864_20250710092726_ch070_000.fits";
    let f1 = args[1].clone();
    let f2 = args[2].clone();
    // wat(f1).unwrap();
    cmp(&f1, &f2)
}
