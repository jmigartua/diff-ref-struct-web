//! Symmetry operators as general (Seitz) operators {R|t}: an integer rotation matrix R acting on fractional
//! coordinates and a translation t, kept exactly in twelfths of a lattice vector and reduced modulo 1.
//!
//! Twelfths are enough for every space-group setting: translations are multiples of 1/6 and Hall's
//! change-of-basis shifts are given in 1/12 (checked over all 530 settings by the test `spglib_530`).
//! A general operator (rotation and translation, not a list of special cases) is what design question Q2
//! asks for, so that magnetic and superspace operators can extend it later.

/// Translations are stored in units of 1/DEN of a lattice vector.
pub const DEN: i32 = 12;

/// A symmetry operator x' = R·x + t with t in units of 1/12, each component in 0..12.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SymOp {
    pub r: [[i32; 3]; 3],
    pub t: [i32; 3],
}

/// The kind of a rotation (proper or improper), from its determinant and trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RotType {
    One,
    Inversion,
    Two,
    Mirror,
    Three,
    InvThree,
    Four,
    InvFour,
    Six,
    InvSix,
}

impl RotType {
    /// The usual symbol: 1, -1, 2, m, 3, -3, 4, -4, 6, -6.
    pub fn symbol(self) -> &'static str {
        match self {
            RotType::One => "1",
            RotType::Inversion => "-1",
            RotType::Two => "2",
            RotType::Mirror => "m",
            RotType::Three => "3",
            RotType::InvThree => "-3",
            RotType::Four => "4",
            RotType::InvFour => "-4",
            RotType::Six => "6",
            RotType::InvSix => "-6",
        }
    }
}

pub const IDENTITY: [[i32; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

impl SymOp {
    pub fn new(r: [[i32; 3]; 3], t: [i32; 3]) -> SymOp {
        SymOp {
            r,
            t: [
                t[0].rem_euclid(DEN),
                t[1].rem_euclid(DEN),
                t[2].rem_euclid(DEN),
            ],
        }
    }

    pub fn identity() -> SymOp {
        SymOp::new(IDENTITY, [0, 0, 0])
    }

    /// A pure translation, e.g. a centring vector.
    pub fn translation(t: [i32; 3]) -> SymOp {
        SymOp::new(IDENTITY, t)
    }

    /// The product self·other: first `other`, then `self`. {R1|t1}{R2|t2} = {R1R2 | R1t2 + t1}.
    pub fn compose(&self, other: &SymOp) -> SymOp {
        let mut r = [[0; 3]; 3];
        let mut t = self.t;
        for (i, row) in r.iter_mut().enumerate() {
            for (j, x) in row.iter_mut().enumerate() {
                *x = (0..3).map(|k| self.r[i][k] * other.r[k][j]).sum();
            }
            t[i] += (0..3).map(|j| self.r[i][j] * other.t[j]).sum::<i32>();
        }
        SymOp::new(r, t)
    }

    pub fn det(&self) -> i32 {
        let m = &self.r;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    pub fn trace(&self) -> i32 {
        self.r[0][0] + self.r[1][1] + self.r[2][2]
    }

    /// The kind of rotation; None if R is not a crystallographic rotation.
    pub fn rot_type(&self) -> Option<RotType> {
        use RotType::*;
        Some(match (self.det(), self.trace()) {
            (1, 3) => One,
            (1, -1) => Two,
            (1, 0) => Three,
            (1, 1) => Four,
            (1, 2) => Six,
            (-1, -3) => Inversion,
            (-1, 1) => Mirror,
            (-1, 0) => InvThree,
            (-1, -1) => InvFour,
            (-1, -2) => InvSix,
            _ => return None,
        })
    }

    /// The operator as coordinate triplet, e.g. "-x+1/2,y,-z+1/2" (FullProf's and CIF's form).
    pub fn to_xyz(&self) -> String {
        (0..3)
            .map(|i| {
                let mut s = String::new();
                for (j, v) in ["x", "y", "z"].iter().enumerate() {
                    let c = self.r[i][j];
                    if c == 0 {
                        continue;
                    }
                    if c < 0 {
                        s.push('-');
                    } else if !s.is_empty() {
                        s.push('+');
                    }
                    if c.abs() > 1 {
                        s.push_str(&c.abs().to_string());
                    }
                    s.push_str(v);
                }
                if self.t[i] != 0 {
                    let g = gcd(self.t[i], DEN);
                    if !s.is_empty() {
                        s.push('+');
                    }
                    s.push_str(&format!("{}/{}", self.t[i] / g, DEN / g));
                } else if s.is_empty() {
                    s.push('0');
                }
                s
            })
            .collect::<Vec<_>>()
            .join(",")
    }

    /// Reads a coordinate triplet such as "x+1/2,-y+1/2,-z" or "-x+y, -x, z+1/3" (case and spaces ignored).
    pub fn parse_xyz(s: &str) -> Result<SymOp, String> {
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 3 {
            return Err(format!(
                "'{s}': a symmetry operator has three comma-separated parts"
            ));
        }
        let mut r = [[0; 3]; 3];
        let mut t = [0; 3];
        for (i, p) in parts.iter().enumerate() {
            let p: String = p
                .chars()
                .filter(|c| !c.is_whitespace())
                .collect::<String>()
                .to_lowercase();
            let mut rest = p.as_str();
            while !rest.is_empty() {
                let (sign, body) = match rest.as_bytes()[0] {
                    b'+' => (1, &rest[1..]),
                    b'-' => (-1, &rest[1..]),
                    _ => (1, rest),
                };
                let end = body.find(['+', '-']).unwrap_or(body.len());
                let term = &body[..end];
                rest = &body[end..];
                match term {
                    "x" => r[i][0] += sign,
                    "y" => r[i][1] += sign,
                    "z" => r[i][2] += sign,
                    _ => {
                        let v = fraction(term)
                            .ok_or_else(|| format!("'{s}': cannot read the term '{term}'"))?;
                        t[i] += sign * v;
                    }
                }
            }
        }
        Ok(SymOp::new(r, t))
    }
}

/// A number such as "1/2", "0.5", "0.3333" or "2/3" in twelfths, if it is a multiple of 1/12.
fn fraction(s: &str) -> Option<i32> {
    let v = match s.split_once('/') {
        Some((a, b)) => a.parse::<f64>().ok()? / b.parse::<f64>().ok()?,
        None => s.parse::<f64>().ok()?,
    };
    let n = (v * DEN as f64).round();
    ((v * DEN as f64 - n).abs() < 1e-3).then_some(n as i32)
}

fn gcd(a: i32, b: i32) -> i32 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn product_and_reduction() {
        let a = SymOp::parse_xyz("-x+1/2,-y,z+1/2").unwrap();
        let aa = a.compose(&a);
        assert_eq!(aa, SymOp::identity());
        let b = SymOp::parse_xyz("x+1/2, -y+1/2, -z+1/2").unwrap();
        assert_eq!(b.t, [6, 6, 6]);
    }

    #[test]
    fn xyz_round_trip() {
        for s in [
            "x,y,z",
            "-x+1/2,y+1/2,-z",
            "-y,x-y,z+1/3",
            "x-y,x,z+5/6",
            "-x+y,y,-z+1/2",
            "z,x,y",
            "x+3/4,y+1/4,z+1/4",
        ] {
            let op = SymOp::parse_xyz(s).unwrap();
            assert_eq!(
                SymOp::parse_xyz(&op.to_xyz()).unwrap(),
                op,
                "{s} -> {}",
                op.to_xyz()
            );
        }
        assert_eq!(
            SymOp::parse_xyz("-x+1/2,y+1/2,-z").unwrap().to_xyz(),
            "-x+1/2,y+1/2,-z"
        );
        assert_eq!(
            SymOp::parse_xyz("X, Y-x, 0.5+Z").unwrap().to_xyz(),
            "x,-x+y,z+1/2"
        );
    }

    #[test]
    fn rotation_types() {
        let t = |s: &str| SymOp::parse_xyz(s).unwrap().rot_type().unwrap().symbol();
        assert_eq!(t("x,y,z"), "1");
        assert_eq!(t("-x,-y,-z"), "-1");
        assert_eq!(t("-x,-y,z"), "2");
        assert_eq!(t("x,y,-z"), "m");
        assert_eq!(t("-y,x-y,z"), "3");
        assert_eq!(t("y,-x+y,-z"), "-3");
        assert_eq!(t("-y,x,z"), "4");
        assert_eq!(t("y,-x,-z"), "-4");
        assert_eq!(t("x-y,x,z"), "6");
        assert_eq!(t("-x+y,-x,-z"), "-6");
    }
}
