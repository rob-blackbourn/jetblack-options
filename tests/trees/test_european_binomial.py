"""Tests for Barone-Adesi-Whaley"""

from jetblack_options.trees.european_binomial import price
from jetblack_options.numeric_greeks.with_carry import NumericGreeks

from ..utils import is_close_to


def test_price():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.070810696746085),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.5069148400613138),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.865263891560955),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.909262426399531),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7893184818731922),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.34561126171898),
    ]:
        b = r - q
        value = price(is_call, S, K, T, r, b, v, 200)
        assert is_close_to(value, expected, 1e-12)


def test_delta():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8606397006611921),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10014973849180042),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5403040031344286),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.42048543601840294),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.166152892001592),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.7946365471511285),
    ]:
        b = r - q
        ng = NumericGreeks(lambda is_call, S, K, T, r, b,
                           v: price(is_call, S, K, T, r, b, v, 100))
        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_gamma():
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.3290705182007514e-11),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 1.1102230246251565e-11),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 7.552788665377008),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 7.552788665390331),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 4.440892098500626e-12),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -1.7763568394002505e-11),
    ]:
        b = r - q
        ng = NumericGreeks(lambda is_call, S, K, T, r, b,
                           v: price(is_call, S, K, T, r, b, v, 100))
        numeric = ng.gamma(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.4761276405391808),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.4187804088877445),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.032411057874453),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.206432268736381),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -2.431468809334608),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.3457394162046157),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b,
                           v: price(is_call, S, K, T, r, b, v, 100))
        b = r - q
        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.58914730397931),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.589147303972009),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.703128947537234),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.703128947531905),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.34480813834438),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.34480813833894),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b,
                           v: price(is_call, S, K, T, r, b, v, 100))
        b = r - q
        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 41.58039938615765),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.981073820943538),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 25.084660240114463),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -22.47681296695525),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.165255681769345),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -44.15236484597873),
    ]:
        ng = NumericGreeks(lambda is_call, S, K, T, r, b,
                           v: price(is_call, S, K, T, r, b, v, 100))
        b = r - q
        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
