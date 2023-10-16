"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.trees.jarrow_rudd import greeks
from jetblack_options.numeric_greeks import NumericGreeks

from ..utils import is_close_to



def test_price():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.068143516623877),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5042503476664845),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.873668938664346),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.917669916891453),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7869051773732514),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.343200400607568),
        (False, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.071997737972442),
        (False, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5164002396952895),
        (False, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8737974667835564),
        (False, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.040445420072334),
        (False, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.786906995272888),
        (False, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 10.099966374860383),
    ]:
        b = r - q
        value, delta, gamma, theta = greeks(is_european, is_call, S, K, T, r, b, v, 200)
        assert is_close_to(value, expected, 1e-12)

def test_delta():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8481107286092637),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.11267866167652407),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.552230551297539),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.4085588389879602),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.14973732129990314),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.8110520689856848),
    ]:
        b = r - q
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)

def test_gamma():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 3.552713678800501e-11),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -8.881784197001252e-12),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -1.7763568394002505e-11),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 1.3322676295501878e-11),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, -4.440892098500626e-12),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -8.881784197001252e-11),
    ]:
        b = r - q
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        numeric = ng.gamma(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.007134014811947864),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.004236199198519475),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.011120478528180211),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.006116587019552089),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.006185897555865183),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.0014244589988159362),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13911625905501346),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13911799017244567),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2669196140807364),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2669211878238531),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17690923380770202),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17691080755082123),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.41173347442779296),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.06388321955324783),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2567353319718184),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.21888136200921204),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.08150631082041687),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.4416720525587108),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
