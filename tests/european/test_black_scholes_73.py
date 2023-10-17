"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.black_scholes_73 import (
    price,
    delta,
    gamma,
    theta,
    vega,
    rho,
    vanna,
    charm,
    ivol
)
from jetblack_options.numeric_greeks_without_carry import NumericGreeksWithoutCarry

from ..utils import is_close_to

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
    ng = NumericGreeksWithoutCarry(price)

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

        numeric = ng.delta(is_call, S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-6)

def test_gamma():
    ng = NumericGreeksWithoutCarry(price)
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

        numeric = ng.gamma(is_call, S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-6)


def test_theta():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, -0.02718824422986413),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.0011271641065568971),
        (True, 100, 100, 0.1, 6/12, 0.125, -0.02623765345826792),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.00017657333496068564),
        (True, 100, 110, 0.1, 6/12, 0.125, -0.0169367334782138),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.011730454657424153),
    ]:
        analytic = theta(is_call, S, K, T, r, v) / 365
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.theta(is_call, S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-6)


def test_vega():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.07463119536211778),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.07463119536211778),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.2342213463902533),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.2342213463902533),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.25278236210643035),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.25278236210643035),
    ]:
        analytic = vega(S, K, T, r, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vega(is_call, S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-3)


def test_rho():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.44954096009369676),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.026073752156660236),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.33244883411948123),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.14316587813087578),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.15110640966088296),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.3720697738145098),
    ]:
        analytic = rho(is_call, S, K, T, r, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.rho(is_call, S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-4)


def test_vanna():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.012280022470048304),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.012280022470048304),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.013819059437024944),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.013819059437024944),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.015924538086889685),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.015924538086889685),
    ]:
        b = r - q
        analytic = vanna(S, K, T, r, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vanna(is_call, S, K, T, r, v, dS=0.01)
        assert is_close_to(numeric, analytic, 1e-4)


def test_charm():
    ng = NumericGreeksWithoutCarry(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.044945814894341574),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.044945814894341574),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.20201591126159346),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.20201591126159346),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.6035085054564097),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.6035085054564097),
    ]:
        b = r - q
        analytic = charm(is_call, S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.charm(is_call, S, K, T, r, v, dS=0.01)
        assert is_close_to(numeric, analytic, 1e-5)
