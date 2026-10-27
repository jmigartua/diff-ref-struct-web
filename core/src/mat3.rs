//! 3×3 matrices and 3-vectors as plain arrays, row-major: `m[i][j]` is row i, column j.

pub type Mat3 = [[f64; 3]; 3];
pub type Vec3 = [f64; 3];

pub fn det(m: &Mat3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Inverse by cofactors; `None` when the matrix is singular to working precision.
pub fn inverse(m: &Mat3) -> Option<Mat3> {
    let d = det(m);
    let scale = m.iter().flatten().fold(0.0_f64, |s, x| s.max(x.abs()));
    if d.abs() <= 1e-14 * scale.powi(3) {
        return None;
    }
    let c =
        |r0: usize, r1: usize, c0: usize, c1: usize| m[r0][c0] * m[r1][c1] - m[r0][c1] * m[r1][c0];
    Some([
        [c(1, 2, 1, 2) / d, -c(0, 2, 1, 2) / d, c(0, 1, 1, 2) / d],
        [-c(1, 2, 0, 2) / d, c(0, 2, 0, 2) / d, -c(0, 1, 0, 2) / d],
        [c(1, 2, 0, 1) / d, -c(0, 2, 0, 1) / d, c(0, 1, 0, 1) / d],
    ])
}

pub fn mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    r
}

pub fn mul_vec(m: &Mat3, v: &Vec3) -> Vec3 {
    [0, 1, 2].map(|i| (0..3).map(|k| m[i][k] * v[k]).sum())
}

/// vᵀ M v
pub fn quad(m: &Mat3, v: &Vec3) -> f64 {
    let mv = mul_vec(m, v);
    (0..3).map(|i| v[i] * mv[i]).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_times_matrix_is_identity() {
        let m = [
            [72.1272, 0.0, -4.8828],
            [0.0, 225.8044, 0.0],
            [-4.8828, 0.0, 130.6604],
        ];
        let p = mul(&m, &inverse(&m).unwrap());
        for i in 0..3 {
            for j in 0..3 {
                let e = if i == j { 1.0 } else { 0.0 };
                assert!((p[i][j] - e).abs() < 1e-12, "{p:?}");
            }
        }
    }

    #[test]
    fn singular_matrix_has_no_inverse() {
        assert!(inverse(&[[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [0.0, 0.0, 1.0]]).is_none());
    }
}
