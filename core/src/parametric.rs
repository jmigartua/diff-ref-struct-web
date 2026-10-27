//! Parametric fits of a sequential refinement's results (P54, L317): a refined value (or a derived one, a cell
//! volume) against the series' external parameter, fitted to an equation the user writes, as GSAS-II's parametric
//! fits of sequential results do (its "Sequential results" window, a user equation of the pseudo-variables and the
//! external parameter, fitted by least squares). The equation is plain text: numbers, the variable (T by default),
//! the coefficients to fit (any other name), + − × / ^, parentheses and the functions exp, ln, log10, sqrt, sin, cos,
//! tan, abs, sinh, cosh, tanh; `pi` is π. Examples: `V0*(1 + a*(T - 150) + b*(T - 150)^2)`; a Birch–Murnaghan
//! equation of state in P.
//!
//! The fit is weighted least squares (weights 1/su²), solved by Levenberg–Marquardt with the derivatives of the parsed
//! equation taken exactly (app request 42); the coefficients' su's are √(C_kk χ²_ν), C the inverse normal matrix and χ²_ν = Σ w r²/(n − p), as
//! the refinement's esd's and SciPy's `curve_fit` with `absolute_sigma=False` (the reference of the test
//! `parametric_fits`).

/// A parsed expression.
#[derive(Clone, Debug, PartialEq)]
enum Expr {
    Num(f64),
    /// The variable.
    Var,
    /// A coefficient, by index.
    Coef(usize),
    Neg(Box<Expr>),
    Bin(char, Box<Expr>, Box<Expr>),
    Fun(&'static str, Box<Expr>),
}

const FUNCTIONS: [&str; 12] = ["exp", "ln", "log10", "sqrt", "sin", "cos", "tan", "abs", "sinh", "cosh", "tanh", "log"];

/// An equation of one variable with named coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct Equation {
    pub text: String,
    pub variable: String,
    /// The coefficients, in the order they first appear.
    pub coefficients: Vec<String>,
    expr: Expr,
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
    variable: &'a str,
    coefficients: Vec<String>,
}

impl Parser<'_> {
    fn skip(&mut self) {
        while self.i < self.s.len() && (self.s[self.i] as char).is_whitespace() {
            self.i += 1;
        }
    }
    fn peek(&mut self) -> Option<char> {
        self.skip();
        self.s.get(self.i).map(|&b| b as char)
    }
    fn err(&self, what: &str) -> String {
        format!("{what} at character {}", self.i + 1)
    }
    fn expr(&mut self) -> Result<Expr, String> {
        let mut e = self.term()?;
        while let Some(c @ ('+' | '-')) = self.peek() {
            self.i += 1;
            e = Expr::Bin(c, Box::new(e), Box::new(self.term()?));
        }
        Ok(e)
    }
    fn term(&mut self) -> Result<Expr, String> {
        let mut e = self.unary()?;
        while let Some(c @ ('*' | '/')) = self.peek() {
            self.i += 1;
            e = Expr::Bin(c, Box::new(e), Box::new(self.unary()?));
        }
        Ok(e)
    }
    fn unary(&mut self) -> Result<Expr, String> {
        match self.peek() {
            Some('-') => {
                self.i += 1;
                Ok(Expr::Neg(Box::new(self.unary()?)))
            }
            Some('+') => {
                self.i += 1;
                self.unary()
            }
            _ => self.power(),
        }
    }
    fn power(&mut self) -> Result<Expr, String> {
        let base = self.primary()?;
        if self.peek() == Some('^') {
            self.i += 1;
            // Right-associative, and binding tighter than a unary minus on its left: -T^2 is -(T^2).
            let exp = self.unary()?;
            return Ok(Expr::Bin('^', Box::new(base), Box::new(exp)));
        }
        Ok(base)
    }
    fn primary(&mut self) -> Result<Expr, String> {
        match self.peek() {
            Some('(') => {
                self.i += 1;
                let e = self.expr()?;
                if self.peek() != Some(')') {
                    return Err(self.err("')' expected"));
                }
                self.i += 1;
                Ok(e)
            }
            Some(c) if c.is_ascii_digit() || c == '.' => {
                let start = self.i;
                while self.i < self.s.len() && ((self.s[self.i] as char).is_ascii_digit() || self.s[self.i] == b'.') {
                    self.i += 1;
                }
                // An exponent: 1e-5, 2.5E3.
                if self.i < self.s.len() && (self.s[self.i] == b'e' || self.s[self.i] == b'E') {
                    let save = self.i;
                    self.i += 1;
                    if self.i < self.s.len() && (self.s[self.i] == b'+' || self.s[self.i] == b'-') {
                        self.i += 1;
                    }
                    if self.i < self.s.len() && (self.s[self.i] as char).is_ascii_digit() {
                        while self.i < self.s.len() && (self.s[self.i] as char).is_ascii_digit() {
                            self.i += 1;
                        }
                    } else {
                        self.i = save;
                    }
                }
                let t = std::str::from_utf8(&self.s[start..self.i]).unwrap();
                t.parse().map(Expr::Num).map_err(|_| self.err(&format!("'{t}' is not a number")))
            }
            Some(c) if c.is_alphabetic() || c == '_' => {
                let start = self.i;
                while self.i < self.s.len() && ((self.s[self.i] as char).is_alphanumeric() || self.s[self.i] == b'_') {
                    self.i += 1;
                }
                let name = std::str::from_utf8(&self.s[start..self.i]).unwrap().to_string();
                if self.peek() == Some('(') {
                    let f = FUNCTIONS.iter().find(|f| **f == name).ok_or(self.err(&format!("'{name}' is not a function")))?;
                    self.i += 1;
                    let arg = self.expr()?;
                    if self.peek() != Some(')') {
                        return Err(self.err("')' expected"));
                    }
                    self.i += 1;
                    return Ok(Expr::Fun(f, Box::new(arg)));
                }
                if name == self.variable {
                    return Ok(Expr::Var);
                }
                if name == "pi" {
                    return Ok(Expr::Num(std::f64::consts::PI));
                }
                let k = self.coefficients.iter().position(|c| *c == name).unwrap_or_else(|| {
                    self.coefficients.push(name.clone());
                    self.coefficients.len() - 1
                });
                Ok(Expr::Coef(k))
            }
            Some(c) => Err(self.err(&format!("'{c}' unexpected"))),
            None => Err(self.err("the equation ends too soon")),
        }
    }
}

fn eval(e: &Expr, x: f64, c: &[f64]) -> f64 {
    match e {
        Expr::Num(v) => *v,
        Expr::Var => x,
        Expr::Coef(k) => c[*k],
        Expr::Neg(a) => -eval(a, x, c),
        Expr::Bin(op, a, b) => {
            let (a, b) = (eval(a, x, c), eval(b, x, c));
            match op {
                '+' => a + b,
                '-' => a - b,
                '*' => a * b,
                '/' => a / b,
                _ => a.powf(b),
            }
        }
        Expr::Fun(f, a) => {
            let a = eval(a, x, c);
            match *f {
                "exp" => a.exp(),
                "ln" | "log" => a.ln(),
                "log10" => a.log10(),
                "sqrt" => a.sqrt(),
                "sin" => a.sin(),
                "cos" => a.cos(),
                "tan" => a.tan(),
                "abs" => a.abs(),
                "sinh" => a.sinh(),
                "cosh" => a.cosh(),
                _ => a.tanh(),
            }
        }
    }
}

/// The value of `e` and its derivative with respect to coefficient k, exactly (forward mode, the chain rule through
/// each operation and function; app request 42).
fn eval_d(e: &Expr, x: f64, c: &[f64], k: usize) -> (f64, f64) {
    match e {
        Expr::Num(v) => (*v, 0.0),
        Expr::Var => (x, 0.0),
        Expr::Coef(j) => (c[*j], if *j == k { 1.0 } else { 0.0 }),
        Expr::Neg(a) => {
            let (v, d) = eval_d(a, x, c, k);
            (-v, -d)
        }
        Expr::Bin(op, a, b) => {
            let ((u, du), (v, dv)) = (eval_d(a, x, c, k), eval_d(b, x, c, k));
            match op {
                '+' => (u + v, du + dv),
                '-' => (u - v, du - dv),
                '*' => (u * v, du * v + u * dv),
                '/' => (u / v, (du * v - u * dv) / (v * v)),
                _ => {
                    let w = u.powf(v);
                    // d(u^v) = v u^(v−1) du + u^v ln u dv; the second term only where v moves (u ≤ 0 otherwise).
                    let dw = if du != 0.0 { v * u.powf(v - 1.0) * du } else { 0.0 } + if dv != 0.0 { w * u.ln() * dv } else { 0.0 };
                    (w, dw)
                }
            }
        }
        Expr::Fun(f, a) => {
            let (u, du) = eval_d(a, x, c, k);
            let (v, dv) = match *f {
                "exp" => (u.exp(), u.exp()),
                "ln" | "log" => (u.ln(), 1.0 / u),
                "log10" => (u.log10(), 1.0 / (u * std::f64::consts::LN_10)),
                "sqrt" => (u.sqrt(), 0.5 / u.sqrt()),
                "sin" => (u.sin(), u.cos()),
                "cos" => (u.cos(), -u.sin()),
                "tan" => (u.tan(), 1.0 / u.cos().powi(2)),
                "abs" => (u.abs(), u.signum()),
                "sinh" => (u.sinh(), u.cosh()),
                "cosh" => (u.cosh(), u.sinh()),
                _ => (u.tanh(), 1.0 - u.tanh().powi(2)),
            };
            (v, if du != 0.0 { dv * du } else { 0.0 })
        }
    }
}

impl Equation {
    /// The equation's derivative with respect to coefficient k at `x` (exact: `eval_d`).
    pub fn derivative(&self, x: f64, c: &[f64], k: usize) -> f64 {
        eval_d(&self.expr, x, c, k).1
    }

    /// Parse `text` as an equation of `variable`.
    pub fn parse(text: &str, variable: &str) -> Result<Equation, String> {
        let mut p = Parser { s: text.as_bytes(), i: 0, variable, coefficients: Vec::new() };
        let expr = p.expr()?;
        if p.peek().is_some() {
            return Err(p.err("unexpected text"));
        }
        if p.coefficients.is_empty() {
            return Err("the equation has no coefficient to fit".into());
        }
        Ok(Equation { text: text.to_string(), variable: variable.to_string(), coefficients: p.coefficients, expr })
    }

    /// The equation's value at `x` with the coefficients `c` (in `coefficients`' order).
    pub fn value(&self, x: f64, c: &[f64]) -> f64 {
        eval(&self.expr, x, c)
    }
}

/// A fitted equation: each coefficient's value and su, χ², χ²_ν, the points and the cycles run.
#[derive(Clone, Debug, PartialEq)]
pub struct Fit {
    pub coefficients: Vec<(String, f64, f64)>,
    /// Σ w r² and Σ w r²/(n − p).
    pub chi2: f64,
    pub chi2_reduced: f64,
    pub points: usize,
    pub cycles: usize,
    /// The coefficients' covariance, C χ²_ν.
    pub covariance: Vec<Vec<f64>>,
}

/// Fit `eq` to the points (x, y, su; a su of 0 or none counts as 1) from the starting coefficients `start` (by name;
/// a coefficient not named starts at 1). Levenberg–Marquardt on the weighted sum of squares, up to 200 cycles, to a
/// relative change of the sum below 10⁻¹⁴.
pub fn fit(eq: &Equation, points: &[(f64, f64, f64)], start: &[(String, f64)]) -> Result<Fit, String> {
    let p = eq.coefficients.len();
    let n = points.len();
    if n <= p {
        return Err(format!("{n} points for {p} coefficients"));
    }
    let mut c: Vec<f64> = eq.coefficients.iter().map(|name| start.iter().find(|s| &s.0 == name).map_or(1.0, |s| s.1)).collect();
    let w: Vec<f64> = points.iter().map(|q| if q.2 > 0.0 && q.2.is_finite() { 1.0 / (q.2 * q.2) } else { 1.0 }).collect();
    let sum = |c: &[f64]| -> f64 { points.iter().zip(&w).map(|(q, wi)| wi * (q.1 - eq.value(q.0, c)).powi(2)).sum() };
    // The derivatives exact, from the parsed equation: a central difference with a step scaled by the coefficient
    // is all rounding when the coefficient fits to 0 (app request 42).
    let jac = |c: &[f64]| -> Vec<Vec<f64>> { (0..p).map(|k| points.iter().map(|q| eq.derivative(q.0, c, k)).collect()).collect() };
    let normal = |c: &[f64]| -> (Vec<Vec<f64>>, Vec<f64>) {
        let j = jac(c);
        let mut a = vec![vec![0.0; p]; p];
        let mut g = vec![0.0; p];
        for i in 0..n {
            let r = points[i].1 - eq.value(points[i].0, c);
            for k in 0..p {
                g[k] += w[i] * j[k][i] * r;
                for l in 0..p {
                    a[k][l] += w[i] * j[k][i] * j[l][i];
                }
            }
        }
        (a, g)
    };
    let mut s = sum(&c);
    if !s.is_finite() {
        return Err("the equation is not finite at the starting coefficients".into());
    }
    let mut lambda = 1e-3;
    let mut cycles = 0;
    for _ in 0..200 {
        cycles += 1;
        let (a, g) = normal(&c);
        let mut improved = false;
        for _ in 0..30 {
            let mut m = a.clone();
            for (k, row) in m.iter_mut().enumerate() {
                row[k] *= 1.0 + lambda;
            }
            let Some(delta) = solve(&m, &g) else {
                lambda *= 10.0;
                continue;
            };
            let trial: Vec<f64> = c.iter().zip(&delta).map(|(x, d)| x + d).collect();
            let st = sum(&trial);
            if st.is_finite() && st <= s {
                let change = (s - st) / s.max(1e-300);
                c = trial;
                s = st;
                lambda = (lambda / 10.0).max(1e-12);
                improved = true;
                if change < 1e-14 {
                    return finish(eq, &c, s, n, cycles, &normal);
                }
                break;
            }
            lambda *= 10.0;
        }
        if !improved {
            break;
        }
    }
    finish(eq, &c, s, n, cycles, &normal)
}

/// The normal matrix and gradient at a set of coefficients.
type Normal<'a> = dyn Fn(&[f64]) -> (Vec<Vec<f64>>, Vec<f64>) + 'a;

fn finish(eq: &Equation, c: &[f64], s: f64, n: usize, cycles: usize, normal: &Normal) -> Result<Fit, String> {
    let p = c.len();
    let (a, _) = normal(c);
    let inv = invert(&a).ok_or("the normal matrix is singular: a coefficient the points cannot fix")?;
    let nu = s / (n - p) as f64;
    let covariance: Vec<Vec<f64>> = inv.iter().map(|r| r.iter().map(|v| v * nu).collect()).collect();
    let coefficients = eq.coefficients.iter().enumerate().map(|(k, name)| (name.clone(), c[k], covariance[k][k].max(0.0).sqrt())).collect();
    Ok(Fit { coefficients, chi2: s, chi2_reduced: nu, points: n, cycles, covariance })
}

/// Gauss–Jordan with partial pivoting.
fn invert(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a.iter().enumerate().map(|(i, r)| r.iter().cloned().chain((0..n).map(|j| if i == j { 1.0 } else { 0.0 })).collect()).collect();
    let scale = a.iter().flatten().fold(0.0f64, |x, v| x.max(v.abs()));
    for col in 0..n {
        let piv = (col..n).max_by(|&i, &j| m[i][col].abs().total_cmp(&m[j][col].abs()))?;
        if m[piv][col].abs() <= 1e-15 * scale {
            return None;
        }
        m.swap(col, piv);
        let d = m[col][col];
        for v in m[col].iter_mut() {
            *v /= d;
        }
        for i in 0..n {
            if i != col {
                let f = m[i][col];
                if f != 0.0 {
                    let row = m[col].clone();
                    for (v, r) in m[i].iter_mut().zip(&row) {
                        *v -= f * r;
                    }
                }
            }
        }
    }
    Some(m.into_iter().map(|r| r[n..].to_vec()).collect())
}

fn solve(a: &[Vec<f64>], g: &[f64]) -> Option<Vec<f64>> {
    let inv = invert(a)?;
    Some(inv.iter().map(|r| r.iter().zip(g).map(|(x, y)| x * y).sum()).collect())
}

/// A column of a series CSV (`sequential::series_csv`) against the external parameter: (parameter, value, su) for
/// every row that has both, the su from the column "<name> su" when there is one.
pub fn column(csv: &str, name: &str, x_column: &str) -> Result<Vec<(f64, f64, f64)>, String> {
    let mut lines = csv.lines();
    let header = split(lines.next().ok_or("an empty table")?);
    let at = |n: &str| header.iter().position(|h| h == n);
    let y = at(name).ok_or(format!("no column '{name}'"))?;
    let x = at(x_column).ok_or(format!("no column '{x_column}'"))?;
    let su = at(&format!("{name} su"));
    let mut out = Vec::new();
    for l in lines {
        let f = split(l);
        let num = |k: usize| f.get(k).and_then(|v| v.parse::<f64>().ok());
        if let (Some(xv), Some(yv)) = (num(x), num(y)) {
            out.push((xv, yv, su.and_then(num).unwrap_or(0.0)));
        }
    }
    Ok(out)
}

/// A CSV line's fields, a quoted field's commas and doubled quotes read as text.
fn split(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '"' if quoted && chars.peek() == Some(&'"') => {
                cur.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ',' if !quoted => out.push(std::mem::take(&mut cur)),
            _ => cur.push(ch),
        }
    }
    out.push(cur);
    out
}
