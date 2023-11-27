"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.trees.leisen_reimer import price
from jetblack_options.numeric_greeks.with_carry import NumericGreeks

from ..utils import is_close_to


def test_price():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.06954773644226),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5056518797578835),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8694961456354164),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.913494680474701),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7881689410025371),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.344461720850616),
        (False, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.073269986822545),
        (False, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5180579183933842),
        (False, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8696227115482875),
        (False, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.036361474985137),
        (False, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7881709020408987),
        (False, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 10.098461348219312),
    ]:
        b = r - q
        value = price(is_european, is_call, S, K, T, r, b, v, 200)
        assert is_close_to(value, expected, 1e-12)


def test_delta():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8567418902401869),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10404754889346002),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5404524757234297),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.42033696347956173),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17152735587044887),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.789262083261999),
    ]:
        b = r - q
        ng = NumericGreeks(
            lambda is_call, S, K, T, r, b, v:
                price(is_european, is_call, S, K, T, r, b, v, 100)
        )
        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_gamma():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018373773116309167),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.01837377170410548),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.04283270254479277),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042832698774475375),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376141771868646),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376155452036755),
    ]:
        b = r - q
        ng = NumericGreeks(
            lambda is_call, S, K, T, r, b, v:
                price(is_european, is_call, S, K, T, r, b, v, 100)
        )
        numeric = ng.gamma(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.514789402139761),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.4574421707461127),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.0402022219039875),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.2142234329979513),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, -2.481116330041932),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.29609189543283954),
    ]:
        ng = NumericGreeks(
            lambda is_call, S, K, T, r, b, v:
                price(is_european, is_call, S, K, T, r, b, v, 100)
        )
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.895002260928102),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.895002260859018),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.769842950124854),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.769842949730947),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.734734325637714),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.734734326296575),
    ]:
        ng = NumericGreeks(
            lambda is_call, S, K, T, r, b, v:
                price(is_european, is_call, S, K, T, r, b, v, 100)
        )
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 41.58596236775036),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.9755108390410205),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 25.08785903306365),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -22.473614173770695),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.182336548597013),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -44.13528397923017),
    ]:
        ng = NumericGreeks(
            lambda is_call, S, K, T, r, b, v:
                price(is_european, is_call, S, K, T, r, b, v, 100)
        )
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
