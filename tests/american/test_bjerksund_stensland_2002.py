"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.american.bjerksund_stensland_2002 import (
    price,
    make_numeric_greeks
)

from ..utils import is_close_to

ng = {
    is_call: make_numeric_greeks(is_call)
    for is_call in (True, False)
}


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
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.857387269854204),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10603562406572564),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5404729148441589),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.44432538850074366),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17153035196955102),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.9106399407253107),
    ]:
        b = r - q
        numeric = ng[is_call].delta(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_gamma():
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018553563609913226),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.01894102808819298),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042840733769367034),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.04794256739160119),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376564245036207),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.04379615347716026),
    ]:
        b = r - q
        numeric = ng[is_call].gamma(S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.530993738452243),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.497522525972883),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.040856526609158),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.4950150852994923),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -2.4811531657213237),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.5957097928179778),
    ]:
        b = r - q
        numeric = ng[is_call].theta(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.978277991237853),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 14.07649727730842),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.773915325762232),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.939020685002646),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.734971740807737),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.059857617765843),
    ]:
        b = r - q
        numeric = ng[is_call].vega(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 40.97350888627194),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.692145861011966),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 25.06962040217786),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -18.115786191373218),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.1822713879518),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -10.273446554156962),
    ]:
        b = r - q
        numeric = ng[is_call].rho(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
