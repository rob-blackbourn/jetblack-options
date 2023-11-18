"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.black_76 import (
    price,
    delta,
    gamma,
    theta,
    vega,
    rho,
    vanna,
    vomma,
    ivol
)
from jetblack_options.numeric_greeks.without_carry import NumericGreeks

from ..utils import is_close_to


def test_price():

    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 10.143390791460092),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.6310965464529535),
        (True, 100, 100, 0.1, 6/12, 0.125, 3.3531192847248605),
        (False, 100, 100, 0.1, 6/12, 0.125, 3.3531192847248534),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.6310965464529654),
        (False, 100, 110, 0.1, 6/12, 0.125, 10.143390791460092),
    ]:
        actual = price(is_call, F, K, T, r, v)
        assert is_close_to(actual, expected, 1e-12)


def test_ivol():

    for is_call, F, K, r, T, p, expected in [
        (True, 110, 100, 0.1, 6/12, 10.143390791460092, 0.125),
        (False, 110, 100, 0.1, 6/12, 0.6310965464529535, 0.125),
        (True, 100, 100, 0.1, 6/12, 3.3531192847248605, 0.125),
        (False, 100, 100, 0.1, 6/12, 3.3531192847248534, 0.125),
        (True, 100, 110, 0.1, 6/12, 0.6310965464529654, 0.125),
        (False, 100, 110, 0.1, 6/12, 10.143390791460092, 0.125),
    ]:
        actual = ivol(is_call, F, K, T, r, p)
        assert is_close_to(actual, expected, 1e-9)


def test_delta():
    ng = NumericGreeks(price)

    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.8267860441956819),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.12444338030503208),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.4923803086739813),
        (False, 100, 100, 0.1, 6/12, 0.125, -0.45884911582673277),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.14319868380006504),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.808030740700649),
    ]:
        analytic = delta(is_call, F, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numerical = ng.delta(is_call, F, K, T, r, v)
        assert is_close_to(numerical, analytic, 1e-6)


def test_gamma():
    ng = NumericGreeks(price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 0.020787293603316447),
        (False, 110, 100, 0.1, 6/12, 0.125, 0.020787293603316447),
        (True, 100, 100, 0.1, 6/12, 0.125, 0.042891991459829734),
        (False, 100, 100, 0.1, 6/12, 0.125, 0.042891991459829734),
        (True, 100, 110, 0.1, 6/12, 0.125, 0.025152625260012922),
        (False, 100, 110, 0.1, 6/12, 0.125, 0.025152625260012922),
    ]:
        analytic = gamma(F, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numerical = ng.gamma(is_call, F, K, T, r, v)
        assert is_close_to(numerical, analytic, 1e-6)


def test_theta():
    ng = NumericGreeks(price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, -0.9507097692924997),
        (False, 110, 100, 0.1, 6/12, 0.125, -1.9019391937932124),
        (True, 100, 100, 0.1, 6/12, 0.125, -3.0156249043267125),
        (False, 100, 100, 0.1, 6/12, 0.125, -3.015624904326714),
        (True, 100, 110, 0.1, 6/12, 0.125, -1.9019391937932129),
        (False, 100, 110, 0.1, 6/12, 0.125, -0.9507097692924997),
    ]:
        analytic = theta(is_call, F, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numerical = ng.theta(is_call, F, K, T, r, v)
        assert is_close_to(numerical, analytic, 1e-4)


def test_vega():
    ng = NumericGreeks(price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, 15.720390787508065),
        (False, 110, 100, 0.1, 6/12, 0.125, 15.720390787508065),
        (True, 100, 100, 0.1, 6/12, 0.125, 26.80749466239359),
        (False, 100, 100, 0.1, 6/12, 0.125, 26.80749466239359),
        (True, 100, 110, 0.1, 6/12, 0.125, 15.720390787508076),
        (False, 100, 110, 0.1, 6/12, 0.125, 15.720390787508076),
    ]:
        analytic = vega(F, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numerical = ng.vega(is_call, F, K, T, r, v, dv=0.005, method='central')
        assert is_close_to(numerical, analytic, 1e-2)


def test_rho():
    ng = NumericGreeks(price)
    for is_call, F, K, r, T, v, expected in [
        (True, 110, 100, 0.1, 6/12, 0.125, -5.071695395730046),
        (False, 110, 100, 0.1, 6/12, 0.125, -0.31554827322647677),
        (True, 100, 100, 0.1, 6/12, 0.125, -1.6765596423624303),
        (False, 100, 100, 0.1, 6/12, 0.125, -1.676559642362427),
        (True, 100, 110, 0.1, 6/12, 0.125, -0.3155482732264827),
        (False, 100, 110, 0.1, 6/12, 0.125, -5.071695395730046),
    ]:
        analytic = rho(is_call, F, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numerical = ng.rho(is_call, F, K, T, r, v)
        assert is_close_to(numerical, analytic, 1e-6)


def test_vanna():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.6720354862986977),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.6720354862986977),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.13403747331196794),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.13403747331196794),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 1.996442942803649),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 1.996442942803649),
    ]:
        b = r - q
        analytic = vanna(S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vanna(is_call, S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-4)


def test_vomma():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 145.98618448189166),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 145.98618448189166),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.41886710409989986),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.41886710409989986),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 145.98618448189166),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 145.98618448189166),
    ]:
        b = r - q
        analytic = vomma(S, K, T, r, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vomma(is_call, S, K, T, r, v)
        assert is_close_to(numeric, analytic, 1e-2)
