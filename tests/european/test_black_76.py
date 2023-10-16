"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.black_76 import (
    price,
)
from jetblack_options.numeric_greeks_without_carry import NumericGreeksWithoutCarry

from ..utils import is_close_to

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
    ng = NumericGreeksWithoutCarry(price)

    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.8267860010469086),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.12444342345385717),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.4923802979506453),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.45884912655054233),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.14319872865428684),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.8080306958468952),
    ]:
        actual = ng.delta(is_call, F, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)

def test_gamma():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.020787294783275456),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.020787294776614118),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.042891986837823026),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.04289198684670481),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.025152624695268244),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.02515262494284798),
    ]:
        actual = ng.gamma(is_call, F, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_theta():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, -0.002604574760095346),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.005211039806134132),
        (True, 100, 100, 0.1, 6/12, 0.125, -0.00827701464981967),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.008277014649812564),
        (True, 100, 110, 0.1, 6/12, 0.125, -0.00521103980614257),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.002604574760095346),
    ]:
        actual = ng.theta(is_call, F, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_vega():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.15684488854121792),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.15684488854121406),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.268074389226308),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.268074389226308),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.15684488854121406),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.15684488854121081),
    ]:
        actual = ng.vega(is_call, F, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_rho():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, -0.05071716527820591),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.0031554958801259247),
        (True, 100, 100, 0.1, 6/12, 0.125, -0.016765666280363245),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.016765666280363245),
        (True, 100, 110, 0.1, 6/12, 0.125, -0.0031554958801259803),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.05071716527820591),
    ]:
        actual = ng.rho(is_call, F, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)
