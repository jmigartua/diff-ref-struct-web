//! Forward-mode dual numbers: a value and its N partial derivatives, carried through arithmetic and the elementary
//! functions exactly (no finite differences). The least-squares engine evaluates the peak positions, widths,
//! intensities and profiles with them, so that its derivatives are analytic (L053).
//!
//! `Real` is the scalar the generic calculations are written for: `f64` for the plain calculation, `Dual<N>` for its
//! derivatives with respect to N chosen variables.

use std::ops::{Add, Div, Mul, Neg, Sub};

/// A scalar the calculations can be written over: `f64` or `Dual<N>`.
pub trait Real:
    Copy + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> + Neg<Output = Self> + Add<f64, Output = Self> + Sub<f64, Output = Self> + Mul<f64, Output = Self> + Div<f64, Output = Self>
{
    fn cst(v: f64) -> Self;
    fn val(self) -> f64;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn tan(self) -> Self;
    fn asin(self) -> Self;
    fn acos(self) -> Self;
    fn sqrt(self) -> Self;
    fn exp(self) -> Self;
    fn ln(self) -> Self;
    fn powi(self, n: i32) -> Self;
    fn powf(self, p: f64) -> Self;
    /// A function of self with value f and derivative df (the chain rule, for functions computed outside).
    fn map(self, f: f64, df: f64) -> Self;
    fn abs(self) -> Self {
        if self.val() < 0.0 {
            -self
        } else {
            self
        }
    }
    fn max0(self) -> Self {
        if self.val() < 0.0 {
            Self::cst(0.0)
        } else {
            self
        }
    }
    fn to_radians(self) -> Self {
        self * (std::f64::consts::PI / 180.0)
    }
    fn to_degrees(self) -> Self {
        self * (180.0 / std::f64::consts::PI)
    }
}

impl Real for f64 {
    fn cst(v: f64) -> Self {
        v
    }
    fn val(self) -> f64 {
        self
    }
    fn sin(self) -> Self {
        f64::sin(self)
    }
    fn cos(self) -> Self {
        f64::cos(self)
    }
    fn tan(self) -> Self {
        f64::tan(self)
    }
    fn asin(self) -> Self {
        f64::asin(self)
    }
    fn acos(self) -> Self {
        f64::acos(self)
    }
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }
    fn exp(self) -> Self {
        f64::exp(self)
    }
    fn ln(self) -> Self {
        f64::ln(self)
    }
    fn powi(self, n: i32) -> Self {
        f64::powi(self, n)
    }
    fn powf(self, p: f64) -> Self {
        f64::powf(self, p)
    }
    fn map(self, f: f64, _df: f64) -> Self {
        f
    }
}

/// A value and its derivatives with respect to N variables.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dual<const N: usize> {
    pub v: f64,
    pub d: [f64; N],
}

impl<const N: usize> Dual<N> {
    /// The variable number i, with value v.
    pub fn var(v: f64, i: usize) -> Self {
        let mut d = [0.0; N];
        d[i] = 1.0;
        Dual { v, d }
    }
    /// A function of self with value f and derivative df (chain rule).
    fn chain(self, f: f64, df: f64) -> Self {
        let mut d = self.d;
        for x in d.iter_mut() {
            *x *= df;
        }
        Dual { v: f, d }
    }
}

impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        let mut d = self.d;
        for (a, b) in d.iter_mut().zip(o.d) {
            *a += b;
        }
        Dual { v: self.v + o.v, d }
    }
}
impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        let mut d = self.d;
        for (a, b) in d.iter_mut().zip(o.d) {
            *a -= b;
        }
        Dual { v: self.v - o.v, d }
    }
}
impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let mut d = [0.0; N];
        for (i, x) in d.iter_mut().enumerate() {
            *x = self.d[i] * o.v + self.v * o.d[i];
        }
        Dual { v: self.v * o.v, d }
    }
}
impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let inv = 1.0 / o.v;
        let mut d = [0.0; N];
        for (i, x) in d.iter_mut().enumerate() {
            *x = (self.d[i] - self.v * inv * o.d[i]) * inv;
        }
        Dual { v: self.v * inv, d }
    }
}
impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    fn neg(self) -> Self {
        self.chain(-self.v, -1.0)
    }
}
impl<const N: usize> Add<f64> for Dual<N> {
    type Output = Self;
    fn add(self, o: f64) -> Self {
        Dual { v: self.v + o, d: self.d }
    }
}
impl<const N: usize> Sub<f64> for Dual<N> {
    type Output = Self;
    fn sub(self, o: f64) -> Self {
        Dual { v: self.v - o, d: self.d }
    }
}
impl<const N: usize> Mul<f64> for Dual<N> {
    type Output = Self;
    fn mul(self, o: f64) -> Self {
        self.chain(self.v * o, o)
    }
}
impl<const N: usize> Div<f64> for Dual<N> {
    type Output = Self;
    fn div(self, o: f64) -> Self {
        self.chain(self.v / o, 1.0 / o)
    }
}

impl<const N: usize> Real for Dual<N> {
    fn cst(v: f64) -> Self {
        Dual { v, d: [0.0; N] }
    }
    fn val(self) -> f64 {
        self.v
    }
    fn sin(self) -> Self {
        self.chain(self.v.sin(), self.v.cos())
    }
    fn cos(self) -> Self {
        self.chain(self.v.cos(), -self.v.sin())
    }
    fn tan(self) -> Self {
        let t = self.v.tan();
        self.chain(t, 1.0 + t * t)
    }
    fn asin(self) -> Self {
        self.chain(self.v.asin(), 1.0 / (1.0 - self.v * self.v).sqrt())
    }
    fn acos(self) -> Self {
        self.chain(self.v.acos(), -1.0 / (1.0 - self.v * self.v).sqrt())
    }
    fn sqrt(self) -> Self {
        let s = self.v.sqrt();
        self.chain(s, if s > 0.0 { 0.5 / s } else { 0.0 })
    }
    fn exp(self) -> Self {
        let e = self.v.exp();
        self.chain(e, e)
    }
    fn ln(self) -> Self {
        self.chain(self.v.ln(), 1.0 / self.v)
    }
    fn powi(self, n: i32) -> Self {
        self.chain(self.v.powi(n), if n == 0 { 0.0 } else { n as f64 * self.v.powi(n - 1) })
    }
    fn powf(self, p: f64) -> Self {
        self.chain(self.v.powf(p), if p == 0.0 { 0.0 } else { p * self.v.powf(p - 1.0) })
    }
    fn map(self, f: f64, df: f64) -> Self {
        self.chain(f, df)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derivatives_of_elementary_functions() {
        let f = |x: Dual<1>| (x.sin() * x.exp() / (x * x + 1.0)).sqrt() + x.acos() * x.powf(1.5) - x.tan().ln();
        let g = |x: f64| (x.sin() * x.exp() / (x * x + 1.0)).sqrt() + x.acos() * x.powf(1.5) - x.tan().ln();
        let x = 0.37;
        let (h, d) = (1e-6, f(Dual::var(x, 0)).d[0]);
        assert!((d - (g(x + h) - g(x - h)) / (2.0 * h)).abs() < 1e-8, "{d}");
    }
}
