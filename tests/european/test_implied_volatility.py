"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.implied_volatility import ivol

from ..utils import is_close_to

def test_ivol():

    for is_call, S, K, r, q, T, p, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 11.069546131685598, 0.125),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.505650275001452, 0.125),
        (True, 100, 100, 0.1, 0.08, 6/12, 3.8695002999527546, 0.125),
        (False, 100, 100, 0.1, 0.08, 6/12, 2.913498834791845, 0.125),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.7881685580252977, 0.125),
        (False, 100, 110, 0.1, 0.08, 6/12, 9.344461337871536, 0.125),
    ]:
        b = r - q
        actual = ivol(is_call, S, K, T, r, b, p)
        assert is_close_to(actual, expected, 1e-9)
