"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.american.bjerksund_stensland_1993 import price
from jetblack_options.numeric_greeks.with_carry import NumericGreeks

from ..utils import is_close_to


def test_price():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.070181515952816),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5100272464024442),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8695089570482253),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.999829098372267),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7881686046834773),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 10.08179174018241),
    ]:
        b = r - q
        actual = price(is_call, S, K, T, r, b, v)
        assert is_close_to(actual, expected, 1e-12)


def test_delta():
    ng = NumericGreeks(price)

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8569705161634467),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10549822644847495),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5404561878748382),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.4430050293766641),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17153014578425996),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.9142855581178999),
    ]:
        b = r - q
        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_gamma():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018450076488818468),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018827352619155135),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.04283403384874873),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.04811970086393558),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.02837645780573439),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.04419040397607432),
    ]:
        b = r - q
        numeric = ng.gamma(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.5209157014280237),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.4877905052248508),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.040357224729991),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.478558694294488),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -2.481145248342571),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.5718287877996886),
    ]:
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.922782099005104),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 14.022558056630885),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.77076342449425),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.87468477418875),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.734918113994524),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.7106813587727327),
    ]:
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 41.357743785999546),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.748840884479023),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 25.083852980962718),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -18.346434662738886),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.182448949924037),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -10.018685692319451),
    ]:
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
