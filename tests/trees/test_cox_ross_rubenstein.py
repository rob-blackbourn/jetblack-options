"""Tests for Cox-Ross-Rubenstein"""

from jetblack_options.trees.cox_ross_rubinstein import (
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
        value = price(is_european, is_call, S, K, T, r, b, v, 200)
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
        numeric = ng[is_european][is_call].delta(S, K, T, r, b, v)
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
        numeric = ng[is_european][is_call].gamma(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_theta():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.4761276405278343),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.4187804088872178),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.032411057869996),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.2064322687333013),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.34573941621434123),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].theta(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_vega():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.589147303977533),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.589147303972009),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.703128947528132),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.70312894752369),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.34480813834477),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.34480813833894),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].vega(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)


def test_rho():

    for is_european, is_call, S, K, r, q, T, v, expected in [
        (True, True, 110, 100, 0.1, 0.08, 6/12, 0.125, 41.58039938613545),
        (True, False, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.981073820944649),
        (True, True, 100, 100, 0.1, 0.08, 6/12, 0.125, 25.084660240106025),
        (True, False, 100, 100, 0.1, 0.08, 6/12, 0.125, -22.476812966961024),
        (True, True, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.165255681767846),
        (True, False, 100, 110, 0.1, 0.08, 6/12, 0.125, -44.152364845996495),
    ]:
        b = r - q
        numeric = ng[is_european][is_call].rho(S, K, T, r, b, v)
        assert is_close_to(numeric, expected, 1e-12)
