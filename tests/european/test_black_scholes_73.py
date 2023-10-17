"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.black_scholes_73 import (
    price,
    ivol
)
from jetblack_options.numeric_greeks_without_carry import NumericGreeksWithoutCarry

from ..utils import is_close_to

def test_price():

    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 15.066208620179964),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.18915107025137257),
        (True, 100, 100, 0.1, 6/12, 0.125, 6.413154785988965),
        (False, 100, 100, 0.1, 6/12, 0.125, 1.536097236060364),
        (True, 100, 110, 0.1, 6/12, 0.125, 1.7525027662779316),
        (False, 100, 110, 0.1, 6/12, 0.125, 6.387739461356475),
    ]:
        actual = price(is_call, S, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)

def test_ivol():

    for is_call, F, K, r, T, p, expected in [
        (True, 110, 100, 0.1, 6/12, 15.066208620179964, 0.125),
        (False, 110, 100, 0.1, 6/12, 0.18915107025137257, 0.125),
        (True, 100, 100, 0.1, 6/12, 6.413154785988965, 0.125),
        (False, 100, 100, 0.1, 6/12, 1.536097236060364, 0.125),
        (True, 100, 110, 0.1, 6/12, 1.7525027662779316, 0.125),
        (False, 100, 110, 0.1, 6/12, 6.387739461356475, 0.125),
    ]:
        actual = ivol(is_call, F, K, T, r, p)
        assert is_close_to(actual, expected, 1e-6)

def test_delta():
    ng = NumericGreeksWithoutCarry(price)

    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.9543127030283927),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.045687296972740654),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.7290291667565896),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.270970833243922),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.3197378759681513),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.6802621240332485),
    ]:
        actual = ng.delta(is_call, S, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)

def test_gamma():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.0098685899274642),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.0098685899274642),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.03747541356347028),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.037475413776633104),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.04044517414314441),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.04044517410761728),
    ]:
        actual = ng.gamma(is_call, S, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_theta():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, -0.02719455697612716),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.0011299065158141275),
        (True, 100, 100, 0.1, 6/12, 0.125, -0.026248771578821106),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.0001841211184903102),
        (True, 100, 110, 0.1, 6/12, 0.125, -0.016918823610495792),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.011752291895859912),
    ]:
        actual = ng.theta(is_call, S, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_vega():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.07456837095173796),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.07456837095173574),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.23400531601991048),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.23400531601991936),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.25258666537197705),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.2525866653719788),
    ]:
        actual = ng.vega(is_call, S, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_rho():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.4494928313661575),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.026123862614649163),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.332388564388729435),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.14322812959207454),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.151143250145358),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.37203511323352245),
    ]:
        actual = ng.rho(is_call, S, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)
