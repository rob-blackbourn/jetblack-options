"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.black_scholes_73 import (
    price,
    make_bumper,
    ivol,
    delta,
    gamma,
    theta,
    vega,
    rho,
    vanna,
    charm,
    vomma
)
from jetblack_options.numeric_greeks.without_carry import NumericGreeks

from ..utils import is_close_to

ng = {
    is_call: make_bumper(is_call)
    for is_call in (True, False)
}


def test_price():

    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 15.066208620179964),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.18915107025137257),
        (True, 100, 100, 0.1, 6/12, 0.125, 6.413154785988965),
        (False, 100, 100, 0.1, 6/12, 0.125, 1.536097236060364),
        (True, 100, 110, 0.1, 6/12, 0.125, 1.7525027662779316),
        (False, 100, 110, 0.1, 6/12, 0.125, 6.387739461356475),
    ]:
        actual = price(is_call, S, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_ivol():

    for is_call, F, K, r, T, p, expected in [
        (True, 110, 100, 0.1, 6/12, 15.066208620179964, 0.125),
        (False, 110, 100, 0.1, 6/12, 0.18915107025137257, 0.125),
        (True, 100, 100, 0.1, 6/12, 6.413154785988965, 0.125),
        (False, 100, 100, 0.1, 6/12, 1.536097236060364, 0.125),
        (True, 100, 110, 0.1, 6/12, 1.7525027662779316, 0.125),
        (False, 100, 110, 0.1, 6/12, 6.387739461356475, 0.125),
    ]:
        actual = ivol(is_call, F, K, T, r, p)
        assert is_close_to(actual, expected, 1e-6)


def test_delta():

    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.9543127330810848),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.045687266918915226),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.7290292160988521),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.2709707839011479),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.3197378469845452),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.6802621530154548),
    ]:
        analytic = delta(is_call, S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng[is_call].delta(S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-6)


def test_gamma():

    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.009868587816478383),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.009868587816478383),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.037475415422440525),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.037475415422440525),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.040445177937028856),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.040445177937028856),
    ]:
        analytic = gamma(S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng[is_call].gamma(S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-6)


def test_theta():

    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, -9.923709143900409),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.41141489889326743),
        (True, 100, 100, 0.1, 6/12, 0.125, -9.576743512267791),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.06444926726065026),
        (True, 100, 110, 0.1, 6/12, 0.125, -6.181907719548038),
        (False, 100, 110, 0.1, 6/12, 0.125, 4.281615949959816),
    ]:
        analytic = theta(is_call, S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng[is_call].theta(S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-4)


def test_vega():

    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 7.463119536211778),
        (False, 110, 100, 0.1, 6/12, 0.125, 7.463119536211778),
        (True, 100, 100, 0.1, 6/12, 0.125, 23.42213463902533),
        (False, 100, 100, 0.1, 6/12, 0.125, 23.42213463902533),
        (True, 100, 110, 0.1, 6/12, 0.125, 25.278236210643035),
        (False, 100, 110, 0.1, 6/12, 0.125, 25.278236210643035),
    ]:
        analytic = vega(S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng[is_call].vega(S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-3)


def test_rho():

    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 44.954096009369676),
        (False, 110, 100, 0.1, 6/12, 0.125, -2.6073752156660236),
        (True, 100, 100, 0.1, 6/12, 0.125, 33.244883411948123),
        (False, 100, 100, 0.1, 6/12, 0.125, -14.316587813087578),
        (True, 100, 110, 0.1, 6/12, 0.125, 15.110640966088296),
        (False, 100, 110, 0.1, 6/12, 0.125, -37.20697738145098),
    ]:
        analytic = rho(is_call, S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng[is_call].rho(S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-4)


def test_vanna():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.2280022470048304),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.2280022470048304),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -1.3819059437024943),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -1.3819059437024943),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 1.5924538086889684),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 1.5924538086889684),
    ]:
        analytic = vanna(S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng[is_call].vanna(S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-4)


def test_charm():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.044945814894341574),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.044945814894341574),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.20201591126159346),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.20201591126159346),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.6035085054564097),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.6035085054564097),
    ]:
        analytic = charm(is_call, S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng[is_call].charm(S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-5)


def test_vomma():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 161.24953775897959),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 161.24953775897959),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 59.59469382217008),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 59.59469382217008),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 52.74707656927028),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 52.74707656927028),
    ]:
        analytic = vomma(S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng[is_call].vomma(S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-2)
