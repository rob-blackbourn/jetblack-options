"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.american.bjerksund_stensland_1993.analytic import price
from jetblack_options.numeric_greeks import NumericGreeks

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
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5404561878737724),
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
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018450076613163446),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.01882735247704659),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.04283403406191155),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.048119701219206945),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376457663625843),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.04419040331882229),
    ]:
        b = r - q
        numeric = ng.gamma(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.006908952112862465),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.004077124553177214),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.011084732983377421),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.006803575363484526),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.006795883337275654),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.0015630348949535744),
    ]:
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13889948587824996),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13989321496247697),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2676966746013001),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2687277038371043),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17698700675118317),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.08525600852293547),
    ]:
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.40721671714833185),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.05756421677021706),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.25055520939913123),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.1836803308496897),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.08187579888723207),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.10087563154881618),
    ]:
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
