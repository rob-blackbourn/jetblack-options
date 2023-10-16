"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.trees.leisen_reimer import greeks
from jetblack_options.numeric_greeks import NumericGreeks

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
        value, delta, gamma, theta = greeks(is_european, is_call, S, K, T, r, b, v, 200)
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
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
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
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        numeric = ng.gamma(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.00689232406346818),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.003994450870044797),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.01108431518136177),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.006080371327426004),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.00679580422672621),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.0008146046747778968),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13862062510248485),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13862062510214607),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2676869855737867),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2676869855738131),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17698516169681633),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17698516169752931),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.41579266366069856),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.05982403031976413),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.25085507412762964),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.22476161985342435),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.08187881464574048),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.441299548733209),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
