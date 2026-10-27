"""The expected values of the order statistics of n standard normal values, by mpmath's quadrature at 30 digits:
E[X_(i)] = ∫ x n!/((i−1)!(n−i)!) Φ(x)^(i−1) (1 − Φ(x))^(n−i) φ(x) dx. The reference for the test
`refinement_diagnostics` (L272), independent of the core's trapezoidal rule. Usage:
python3 expected_normal_order_statistics.py 10          (the upper half of n = 10)
python3 expected_normal_order_statistics.py 13000 1 101  (orders 1 and 101 of n = 13000)
"""
import sys
import mpmath as mp

mp.mp.dps = 30
n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
orders = [int(a) for a in sys.argv[2:]] or range(n // 2 + 1, n + 1)
for i in orders:
    c = mp.exp(mp.loggamma(n + 1) - mp.loggamma(i) - mp.loggamma(n - i + 1))
    f = lambda x: x * c * mp.ncdf(x) ** (i - 1) * (1 - mp.ncdf(x)) ** (n - i) * mp.npdf(x)
    print(i, mp.nstr(mp.quad(f, [-mp.inf, -6, -4.5, -4, -3.5, -3, -2.5, -2, 0, 2, 3, mp.inf]), 15))
