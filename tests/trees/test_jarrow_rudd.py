"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.trees.jarrow_rudd import (
    price,
    make_bumper
)

from ..utils import is_close_to

ng = {
    is_european: {
        is_call: make_bumper(is_european, is_call, 100)
        for is_call in (True, False)
    }
    for is_european in (True, False)
}


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
        value = price(is_european, is_call, S, K, T, r, b, v, 200)
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
        numeric = ng[is_european][is_call].delta(S, K, T, r, b, v)
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
        numeric = ng[is_european][is_call].gamma(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.600415301930754),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.5430891416259644),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.053377000500915),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.227417367121135),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, -2.254806889349802),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.5223821804422535),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].theta(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 14.732851803984204),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 14.733023825931806),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.69217653459088),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.692332918180515),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 16.268533431813616),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 16.26868981540408),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].vega(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 41.11087648857925),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -6.450596718397261),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 25.6734272942718),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -21.888045912704655),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 7.095656390242011),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -45.221964137427406),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].rho(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
