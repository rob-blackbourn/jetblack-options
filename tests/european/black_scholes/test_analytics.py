"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.black_scholes.analytic import (
    price,
    delta,
    gamma,
    theta,
    vega,
    rho,
)
from jetblack_options.numeric_greeks import NumericGreeks

from ...utils import is_close_to

def test_price():
    S = 110
    K = 100
    r = 0.1
    q = 0.08
    T = 6 / 12
    v = 0.125
    b = r - q
    actual = price(True, S, K, T, r, b, v)
    assert is_close_to(actual, 10.563893298711719, 1e-12)

def test_delta():
    S = 110
    K = 100
    r = 0.1
    q = 0.08
    T = 6 / 12
    v = 0.125
    b = r - q
    actual = delta(True, S, K, T, r, b, v)
    assert is_close_to(actual, 0.9607868160720119, 1e-12)

def test_numeric_delta():
    ng = NumericGreeks(price)
    S = 110
    K = 100
    r = 0.1
    q = 0.08
    T = 6 / 12
    v = 0.125
    b = r - q

    actual = ng.delta(True, S, K, T, r, b, v)
    expected = delta(True, S, K, T, r, b, v)
    assert is_close_to(actual, expected, 1e-5)

def test_numeric_gamma():
    ng = NumericGreeks(price)
    S = 110
    K = 100
    r = 0.1
    q = 0.08
    T = 6 / 12
    v = 0.125
    b = r - q

    actual = ng.gamma(True, S, K, T, r, b, v, dS=0.01)
    expected = gamma(S, K, T, r, b, v)
    assert is_close_to(actual, expected, 1e-5)


def test_numeric_theta():
    ng = NumericGreeks(price)
    is_call = True
    S = 110
    K = 100
    r = 0.1
    q = 0.08
    T = 6 / 12
    v = 0.125
    b = r - q

    actual = ng.theta(is_call, S, K, T, r, b, v)
    expected = theta(is_call, S, K, T, r, b, v) / 365

    assert is_close_to(actual, expected, 1e-5)


def test_numeric_vega():
    ng = NumericGreeks(price)
    is_call = True
    S = 110
    K = 100
    r = 0.1
    q = 0.08
    T = 6 / 12
    v = 0.125
    b = r - q

    actual = ng.vega(is_call, S, K, T, r, b, v)
    expected = vega(S, K, T, r, b, v) / 100

    assert is_close_to(actual, expected, 1e-3)


def test_numeric_rho():
    ng = NumericGreeks(price)
    is_call = True
    S = 110
    K = 100
    r = 0.1
    q = 0.08
    T = 6 / 12
    v = 0.125
    b = r - q

    actual = ng.rho(is_call, S, K, T, r, b, v)
    expected = rho(is_call, S, K, T, r, b, v) / 100

    assert is_close_to(actual, expected, 1e-4)
