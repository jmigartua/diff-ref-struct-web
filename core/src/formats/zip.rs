//! Zip archives, enough for the data formats that are zips (Bruker `.brml`, Rigaku `.rasx`, MIDAS `.zarr.zip`): the
//! central directory (PKWARE APPNOTE 6.3.10, §4.3.12–4.3.16), stored and deflated entries, and DEFLATE itself
//! (RFC 1951). No Zip64, no encryption.

use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct ZipError(pub String);

impl fmt::Display for ZipError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(&self.0)
    }
}

fn zerr<T>(m: impl Into<String>) -> Result<T, ZipError> {
    Err(ZipError(m.into()))
}

fn u16_at(b: &[u8], i: usize) -> Option<u16> {
    Some(u16::from_le_bytes(b.get(i..i + 2)?.try_into().ok()?))
}

fn u32_at(b: &[u8], i: usize) -> Option<u32> {
    Some(u32::from_le_bytes(b.get(i..i + 4)?.try_into().ok()?))
}

/// An entry of the central directory.
#[derive(Clone, Debug)]
pub struct Entry {
    pub name: String,
    method: u16,
    compressed: usize,
    size: usize,
    offset: usize,
}

/// A zip archive held in memory.
pub struct Archive<'a> {
    bytes: &'a [u8],
    pub entries: Vec<Entry>,
}

impl<'a> Archive<'a> {
    /// The central directory, found through the end-of-central-directory record (signature 0x06054b50) in the last
    /// 64 KiB + 22 bytes.
    pub fn open(bytes: &'a [u8]) -> Result<Archive<'a>, ZipError> {
        if bytes.len() < 22 || u32_at(bytes, 0) != Some(0x0403_4b50) {
            return zerr("not a zip archive");
        }
        let lo = bytes.len().saturating_sub(65_557);
        let eocd = (lo..=bytes.len() - 22).rev().find(|&i| u32_at(bytes, i) == Some(0x0605_4b50)).ok_or(ZipError("zip: no end of central directory".into()))?;
        let n = u16_at(bytes, eocd + 10).unwrap_or(0) as usize;
        let mut at = u32_at(bytes, eocd + 16).unwrap_or(u32::MAX) as usize;
        let mut entries = Vec::with_capacity(n.min(4096));
        for _ in 0..n {
            if u32_at(bytes, at) != Some(0x0201_4b50) {
                return zerr("zip: bad central directory entry");
            }
            let get16 = |k| u16_at(bytes, at + k).ok_or(ZipError("zip: truncated directory".into()));
            let get32 = |k| u32_at(bytes, at + k).ok_or(ZipError("zip: truncated directory".into()));
            let (method, compressed, size, offset) = (get16(10)?, get32(20)? as usize, get32(24)? as usize, get32(42)? as usize);
            let (ln, le, lc) = (get16(28)? as usize, get16(30)? as usize, get16(32)? as usize);
            let name = String::from_utf8_lossy(bytes.get(at + 46..at + 46 + ln).ok_or(ZipError("zip: truncated name".into()))?).into_owned();
            if compressed == u32::MAX as usize || offset == u32::MAX as usize {
                return zerr("zip: Zip64 archives are not read");
            }
            entries.push(Entry { name, method, compressed, size, offset });
            at += 46 + ln + le + lc;
        }
        Ok(Archive { bytes, entries })
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.entries.iter().map(|e| e.name.as_str())
    }

    /// The bytes of the entry with this name.
    pub fn read(&self, name: &str) -> Result<Vec<u8>, ZipError> {
        let e = self.entries.iter().find(|e| e.name == name).ok_or(ZipError(format!("zip: no entry {name}")))?;
        let b = self.bytes;
        if u32_at(b, e.offset) != Some(0x0403_4b50) {
            return zerr(format!("zip: bad local header for {name}"));
        }
        let start = e.offset + 30 + u16_at(b, e.offset + 26).unwrap_or(0) as usize + u16_at(b, e.offset + 28).unwrap_or(0) as usize;
        let data = b.get(start..start + e.compressed).ok_or(ZipError(format!("zip: {name} truncated")))?;
        let out = match e.method {
            0 => data.to_vec(),
            8 => inflate(data)?,
            m => return zerr(format!("zip: {name} uses compression method {m}, not read")),
        };
        if out.len() != e.size {
            return zerr(format!("zip: {name} inflates to {} bytes, the directory says {}", out.len(), e.size));
        }
        Ok(out)
    }
}

/// Bits read least-significant first (RFC 1951 §3.1.1).
struct Bits<'a> {
    b: &'a [u8],
    pos: usize,
    bit: u32,
    nbits: u32,
}

impl Bits<'_> {
    fn need(&mut self, n: u32) -> Result<u32, ZipError> {
        while self.nbits < n {
            let byte = *self.b.get(self.pos).ok_or(ZipError("deflate: data ends early".into()))?;
            self.pos += 1;
            self.bit |= (byte as u32) << self.nbits;
            self.nbits += 8;
        }
        let v = self.bit & ((1u64 << n) - 1) as u32;
        self.bit >>= n;
        self.nbits -= n;
        Ok(v)
    }
}

/// A canonical Huffman code (RFC 1951 §3.2.2): the number of codes of each length and the symbols in code order.
struct Huffman {
    count: [u16; 16],
    symbol: Vec<u16>,
}

impl Huffman {
    fn new(lengths: &[u8]) -> Huffman {
        let mut count = [0u16; 16];
        for &l in lengths {
            count[l as usize] += 1;
        }
        count[0] = 0;
        let mut offs = [0u16; 16];
        for l in 1..16 {
            offs[l] = offs[l - 1] + count[l - 1];
        }
        let mut symbol = vec![0u16; lengths.len()];
        for (s, &l) in lengths.iter().enumerate() {
            if l != 0 {
                symbol[offs[l as usize] as usize] = s as u16;
                offs[l as usize] += 1;
            }
        }
        Huffman { count, symbol }
    }

    fn decode(&self, bits: &mut Bits) -> Result<u16, ZipError> {
        let (mut code, mut first, mut index) = (0i32, 0i32, 0i32);
        for l in 1..16 {
            code |= bits.need(1)? as i32;
            let c = self.count[l] as i32;
            if code - c < first {
                return self.symbol.get((index + (code - first)) as usize).copied().ok_or(ZipError("deflate: bad code".into()));
            }
            index += c;
            first = (first + c) << 1;
            code <<= 1;
        }
        zerr("deflate: bad code")
    }
}

const LEN_BASE: [u16; 29] = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
const LEN_EXTRA: [u8; 29] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
const DIST_BASE: [u16; 30] = [1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577];
const DIST_EXTRA: [u8; 30] = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];

fn codes(bits: &mut Bits, out: &mut Vec<u8>, lit: &Huffman, dist: &Huffman) -> Result<(), ZipError> {
    loop {
        let s = lit.decode(bits)? as usize;
        if s < 256 {
            out.push(s as u8);
        } else if s == 256 {
            return Ok(());
        } else {
            let k = s - 257;
            if k >= 29 {
                return zerr("deflate: bad length symbol");
            }
            let len = LEN_BASE[k] as usize + bits.need(LEN_EXTRA[k] as u32)? as usize;
            let d = dist.decode(bits)? as usize;
            if d >= 30 {
                return zerr("deflate: bad distance symbol");
            }
            let back = DIST_BASE[d] as usize + bits.need(DIST_EXTRA[d] as u32)? as usize;
            if back > out.len() {
                return zerr("deflate: distance before the start");
            }
            let from = out.len() - back;
            for i in 0..len {
                out.push(out[from + i]);
            }
        }
    }
}

/// Raw DEFLATE data (RFC 1951) inflated.
pub fn inflate(data: &[u8]) -> Result<Vec<u8>, ZipError> {
    let mut bits = Bits { b: data, pos: 0, bit: 0, nbits: 0 };
    let mut out = Vec::with_capacity(data.len() * 4);
    loop {
        let last = bits.need(1)?;
        match bits.need(2)? {
            0 => {
                // Stored: to the byte boundary, LEN and NLEN, LEN bytes.
                bits.bit = 0;
                bits.nbits = 0;
                let p = bits.pos;
                let (len, nlen) = (u16_at(data, p).ok_or(ZipError("deflate: truncated".into()))?, u16_at(data, p + 2).ok_or(ZipError("deflate: truncated".into()))?);
                if len != !nlen {
                    return zerr("deflate: stored block length check");
                }
                let block = data.get(p + 4..p + 4 + len as usize).ok_or(ZipError("deflate: truncated".into()))?;
                out.extend_from_slice(block);
                bits.pos = p + 4 + len as usize;
            }
            1 => {
                let mut l = [0u8; 288];
                l[..144].fill(8);
                l[144..256].fill(9);
                l[256..280].fill(7);
                l[280..].fill(8);
                codes(&mut bits, &mut out, &Huffman::new(&l), &Huffman::new(&[5u8; 30]))?;
            }
            2 => {
                let (nlen, ndist, ncode) = (bits.need(5)? as usize + 257, bits.need(5)? as usize + 1, bits.need(4)? as usize + 4);
                const ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
                let mut cl = [0u8; 19];
                for &o in ORDER.iter().take(ncode) {
                    cl[o] = bits.need(3)? as u8;
                }
                let ch = Huffman::new(&cl);
                let mut lengths = Vec::with_capacity(nlen + ndist);
                while lengths.len() < nlen + ndist {
                    let s = ch.decode(&mut bits)?;
                    match s {
                        0..=15 => lengths.push(s as u8),
                        16 => {
                            let prev = *lengths.last().ok_or(ZipError("deflate: repeat with no length".into()))?;
                            for _ in 0..3 + bits.need(2)? {
                                lengths.push(prev);
                            }
                        }
                        17 => lengths.extend(std::iter::repeat_n(0, 3 + bits.need(3)? as usize)),
                        _ => lengths.extend(std::iter::repeat_n(0, 11 + bits.need(7)? as usize)),
                    }
                }
                if lengths.len() > nlen + ndist {
                    return zerr("deflate: too many lengths");
                }
                codes(&mut bits, &mut out, &Huffman::new(&lengths[..nlen]), &Huffman::new(&lengths[nlen..]))?;
            }
            _ => return zerr("deflate: block type 3"),
        }
        if last == 1 {
            return Ok(out);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inflate_fixed_and_stored() {
        // Raw DEFLATE of "hello hello hello hello" by zlib (level 9, wbits −15): fixed Huffman, one back-reference.
        let d = [0xcb, 0x48, 0xcd, 0xc9, 0xc9, 0x57, 0xc8, 0x40, 0x27, 0x01];
        assert_eq!(inflate(&d).unwrap(), b"hello hello hello hello");
        // A stored block: BFINAL=1, BTYPE=00, LEN=3, NLEN=!3, "abc".
        assert_eq!(inflate(&[0x01, 0x03, 0x00, 0xfc, 0xff, b'a', b'b', b'c']).unwrap(), b"abc");
    }
}
