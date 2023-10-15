"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.american.barone_adesi_whaley.analytic import price
from jetblack_options.numeric_greeks import NumericGreeks

from ..utils import is_close_to



def test_price():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.087510335081676),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.510639694796271),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8736244925135566),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.938715732901822),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7892100659783038),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.484654220828427),
    ]:
        b = r - q
        actual = price(is_call, S, K, T, r, b, v)
        assert is_close_to(actual, expected, 1e-12)

def test_delta():
    ng = NumericGreeks(price)

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8592614442108903),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10482043749096559),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.541088570403403),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.42462427819012216),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17169091749034138),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.813090973112196),
    ]:
        b = r - q
        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)

def test_gamma():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018705123903117737),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018500320174696938),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.04292392062232864),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.043603545081261075),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.02839965949075207),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.03266594109874177),
    ]:
        b = r - q
        numeric = ng.gamma(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.007141292450844716),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.004052889472088705),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.011158298535876021),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.006272273350176771),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.006818777514447372),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.00032474140609117796),
    ]:
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.1423167545129651),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13950445036474848),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2689702277089201),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.26851816368330694),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17742200338091113),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.16099953881893558),
    ]:
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.3853122903774979),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.057185434099532584),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.24431887356790716),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.2096309209889331),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.08034077274530887),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.3468338033537526),
    ]:
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
