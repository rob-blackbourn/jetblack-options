"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.trees.trinomial import (
    price,
    make_numeric_greeks
)

from ..utils import is_close_to

ng = {
    is_european: {
        is_call: make_numeric_greeks(is_european, is_call, 100)
        for is_call in (True, False)
    }
    for is_european in (True, False)
}


def test_price():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.073487046391277),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5177978858542516),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8675119431883433),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.0348490093800793),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7880335120168521),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 10.099148293480416),
        (False, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.069621888515023),
        (False, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5057260318285591),
        (False, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.867381427247627),
        (False, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.911379962084615),
        (False, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7880316003855702),
        (False, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.344324380229715),
    ]:
        b = r - q
        value = price(is_european, is_call, S, K, T, r, b, v, 200)
        assert is_close_to(value, expected, 1e-12)


def test_delta():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8546747393094023),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10969294624645909),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5404289451924393),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.4507105497704389),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17600965898771914),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.8991782315475483),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].delta(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_gamma():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 5.898801447301594e-05),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0047883994047648315),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 5.347308742291723),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.3960301873525722),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.00746144338882e-08),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.0033522500508809117),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].gamma(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.576238633940049),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.546128615272118),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.037857952421628),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.522968202260064),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, -2.5246895681611035),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.6630983238472599),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].theta(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 14.313819649525605),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 14.411707802939155),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.746272777726343),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.939561719354188),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 18.087281241793207),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 10.159876255306344),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].vega(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 40.46322386017476),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.514069958272582),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 25.037258609427404),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -17.291090496786452),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.17243182565347),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -11.521926366350854),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].rho(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
