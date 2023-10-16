"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.black_76 import (
    price,
)
from jetblack_options.numeric_greeks import NumericGreeks

from ..utils import is_close_to

def _ng_price(
        is_call: bool,
        F: float,
        K: float,
        r: float,
        _b: float,
        T: float,
        v: float
) -> float:
    return price(is_call, F, K, T, r, v)

def test_price():

    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 10.143390791460092),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.6310965464529535),
        (True, 100, 100, 0.1, 6/12, 0.125, 3.3531192847248605),
        (False, 100, 100, 0.1, 6/12, 0.125, 3.3531192847248534),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.6310965464529654),
        (False, 100, 110, 0.1, 6/12, 0.125, 10.143390791460092),
    ]:
        actual = price(is_call, F, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)

def test_delta():
    ng = NumericGreeks(_ng_price)

    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.9440668072311809),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.007162617270356237),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.48311444341376797),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.4681149810870866),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.00798259369692517),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.9432468308049557),
    ]:
        actual = ng.delta(is_call, F, K, T, r, r, v)
        assert is_close_to(actual, expected, 1e-12)

def test_gamma():
    ng = NumericGreeks(_ng_price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.004546386271897518),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.004546386066835856),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.09598431623958348),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.09598431630619686),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.0055011274911580255),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.005501127571250208),
    ]:
        actual = ng.gamma(is_call, F, K, T, r, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_theta():
    ng = NumericGreeks(_ng_price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.0026093067367956024),
        (False, 110, 100, 0.1, 6/12, 0.125, 2.8416907624020937e-06),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.00041100178721920066),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.00041100178721920066),
        (True, 100, 110, 0.1, 6/12, 0.125, 2.8416907624020937e-06),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.0026093067367956024),
    ]:
        actual = ng.theta(is_call, F, K, T, r, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_vega():
    ng = NumericGreeks(_ng_price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.006995853595110013),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.006995853595114785),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.11998040889873018),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.11998040889873018),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.006995853595117426),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.006995853595117119),
    ]:
        actual = ng.vega(is_call, F, K, T, r, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_rho():
    ng = NumericGreeks(_ng_price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, -0.04331478628344421),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.004246883114637223),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.06759251730340632),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.06759251730340632),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.004246883114637013),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.04331478628344421),
    ]:
        actual = ng.rho(is_call, F, K, T, r, r, v)
        assert is_close_to(actual, expected, 1e-12)
