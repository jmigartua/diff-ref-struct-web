//! The two container formats that GSAS-II reads through Python libraries (`h5py`, `zarr`): recognised here by their
//! signatures, not read. Reading them without a dependency needs an HDF5 decoder (superblock, object headers,
//! B-trees, chunked and deflated datasets) and zarr's Blosc codecs (LZ4, Zstandard); the `hdf5` crate binds the HDF5
//! C library, and pure-Rust readers were not evaluated (L062). No example file of either format was found to test a
//! reader against. Both are left for a later release.
//! - NeXus/HDF5 (`G2pwd_HDF5.py`, MAX IV NXazint1d and NXazint2d): the HDF5 signature `\x89HDF\r\n\x1a\n` at byte 0 (or
//!   at 512, 1024, … for a file with a user block).
//! - MIDAS zarr (`G2pwd_MIDAS.py`): a zip holding zarr metadata (`.zgroup`, `.zarray` or `zarr.json`) and the
//!   `REtaMap` and `OmegaSumFrame` arrays.

use super::{Probe, SIGNATURE};

const HDF5: &[u8; 8] = b"\x89HDF\r\n\x1a\n";

pub fn recognise_hdf5(p: &Probe) -> f64 {
    let mut at = 0usize;
    while at + 8 <= p.bytes.len() && at <= 1 << 20 {
        if &p.bytes[at..at + 8] == HDF5 {
            // An HDF5 file; NeXus azimuthal integration files carry their definition name in the file.
            let nx = p.bytes.windows(9).any(|w| w == b"NXazint1d" || w == b"NXazint2d");
            return if nx { SIGNATURE } else { SIGNATURE - 0.3 };
        }
        at = if at == 0 { 512 } else { at * 2 };
    }
    0.0
}

pub fn recognise_zarr(p: &Probe) -> f64 {
    match p.zip_names() {
        Some(n) => {
            let meta = n.iter().any(|e| e.ends_with(".zgroup") || e.ends_with(".zarray") || e.ends_with("zarr.json"));
            let midas = n.iter().any(|e| e.contains("REtaMap")) && n.iter().any(|e| e.contains("OmegaSumFrame"));
            if meta && midas {
                SIGNATURE
            } else {
                0.0
            }
        }
        None => 0.0,
    }
}
