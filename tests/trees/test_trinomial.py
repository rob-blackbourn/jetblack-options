"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.trees.trinomial import greeks
from jetblack_options.numeric_greeks.with_carry import NumericGreeks

from ..utils import is_close_to



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
        value, delta, gamma, theta = greeks(is_european, is_call, S, K, T, r, b, v, 200)
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
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        numeric = ng.delta(is_call, S, K, T, r, b, v)
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
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        numeric = ng.gamma(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.007068073602292557),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.004243682405050797),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.01107782687525205),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.006925358650654356),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.006925015063197537),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.0018188950369175672),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13876296515088526),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.14029914711009464),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2674544353724515),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.26937512505191696),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.175380456387576),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.09277815754806085),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.3916062832273983),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.05525182662608269),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.24885333736768023),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.17323262467686518),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.0817104036732455),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.10971075485336357),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
