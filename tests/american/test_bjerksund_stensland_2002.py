"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.american.bjerksund_stensland_2002 import price
from jetblack_options.numeric_greeks.with_carry import NumericGreeks

from ..utils import is_close_to



def test_price():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.071608152504766),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5121325794818574),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8695471290130214),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.013317901457917),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7881689776961025),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 10.087984869664105),
    ]:
        b = r - q
        actual = price(is_call, S, K, T, r, b, v)
        assert is_close_to(actual, expected, 1e-12)

def test_delta():
    ng = NumericGreeks(price)

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.857387269854204),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10603562406572564),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5404729148441589),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.44432538850074366),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17153035196955102),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.9106399407253107),
    ]:
        b = r - q
        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)

def test_gamma():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018553563609913226),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.01894102808819298),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042840733769367034),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.04794256739160119),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376564245036207),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.04379615347716026),
    ]:
        b = r - q
        numeric = ng.gamma(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.006936412724918739),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.004103743695836215),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.011086082697737254),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.006848792611478416),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.006795904502169492),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.001628723104140306),
    ]:
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.139454331796788),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.1404303714778692),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2677295300944209),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.26937284966896513),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17698765412465178),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.08874969924391607),
    ]:
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.40084245863918877),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.05700346376843157),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2500455606081218),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.18130969095967941),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.08186189553637746),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.10308634314698306),
    ]:
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
