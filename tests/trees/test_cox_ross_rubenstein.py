"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.trees.cox_ross_rubinstein import greeks
from jetblack_options.numeric_greeks import NumericGreeks

from ..utils import is_close_to



def test_price():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.070810696746085),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5069148400613138),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.865263891560955),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.909262426399531),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7893184818731922),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.34561126171898),
        (False, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.074659452363559),
        (False, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5190977256236272),
        (False, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8653919092508406),
        (False, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.034209281414657),
        (False, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7893202834043941),
        (False, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 10.100588559015748),
    ]:
        b = r - q
        value, delta, gamma, theta = greeks(is_european, is_call, S, K, T, r, b, v, 200)
        assert is_close_to(value, expected, 1e-12)

def test_delta():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8606397006611921),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10014973849180042),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5403040031344286),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.42048543601840294),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.166152892001592),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.7946365471511285),
    ]:
        b = r - q
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)

def test_gamma():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.3290705182007514e-11),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 1.6653345369377348e-11),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 7.552788665350363),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 7.5527886653592454),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 1.0658141036401503e-10),
    ]:
        b = r - q
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        numeric = ng.gamma(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.006793772874507553),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.003895899679320025),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.011062933026733557),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.0060589891710591814),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.006669281589761655),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.0009411273119681596),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13611230457699985),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13611230457699058),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2670198719678558),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2670198719678494),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.1739838829453611),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17398388294535572),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.41573822067700394),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.05987847330393839),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.25082336917322356),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.22479332480772052),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.0817070259188224),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.44147133746021705),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
