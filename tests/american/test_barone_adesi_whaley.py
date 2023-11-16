"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.american.barone_adesi_whaley import price
from jetblack_options.numeric_greeks.with_carry import NumericGreeks

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
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.01870512376100919),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.01850032019357073),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042923920551274364),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.04360354515231535),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028399659561806345),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.03266594122308675),
    ]:
        b = r - q
        numeric = ng.gamma(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.6059936396929517),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.4788290399032804),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.067355373044116),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.284344775164011),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -2.48956425165973),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.11716003329799829),
    ]:
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 14.2656452419061),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.9833977312675),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.89826727142619),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.853425641593763),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.77822568607357),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 16.181880739890353),
    ]:
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 39.274542330740125),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.709938455815355),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 24.580558678367613),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -20.964770823997945),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.060271439361),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -34.82234069684775),
    ]:
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
