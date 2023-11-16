"""Tests for generalised Black-Scholes European options"""

from jetblack_options.european.generalised_black_scholes import (
    price,
    ivol,
    delta,
    gamma,
    theta,
    vega,
    rho,
    elasticity,
    dgamma_dvol,
    gammap,
    vanna,
    vegap,
    charm,
    vomma,
)
from jetblack_options.numeric_greeks.with_carry import NumericGreeks

from ..utils import is_close_to


def test_price():

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.069546131685598),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.505650275001452),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8695002999527546),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.913498834791845),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7881685580252977),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.344461337871536),
    ]:
        b = r - q
        actual = price(is_call, S, K, T, r, b, v)
        assert is_close_to(actual, expected, 1e-12)


def test_ivol():

    for is_call, S, K, r, q, T, p, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 11.069546131685598, 0.125),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.505650275001452, 0.125),
        (True, 100, 100, 0.1, 0.08, 6/12, 3.8695002999527546, 0.125),
        (False, 100, 100, 0.1, 0.08, 6/12, 2.913498834791845, 0.125),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.7881685580252977, 0.125),
        (False, 100, 110, 0.1, 0.08, 6/12, 9.344461337871536, 0.125),
    ]:
        b = r - q
        actual = ivol(is_call, S, K, T, r, b, p)
        assert is_close_to(actual, expected, 1e-9)


def test_delta():
    ng = NumericGreeks(price)

    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8567400985874144),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10404934056490878),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5404518486173583),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.42033759053496483),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17153007262292186),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.7892593665294013),
    ]:
        b = r - q
        analytic = delta(is_call, S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, 1e-5)


def test_gamma():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018374151835767315),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018374151835767315),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042831984686328525),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042831984686328525),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376442324910798),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376442324910798),
    ]:
        b = r - q
        analytic = gamma(S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.gamma(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, analytic, 1e-5)


def test_theta():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -2.514805144448628),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -1.457457963981934),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -4.040202470858858),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -2.2142237390703023),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -2.48115284607697),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.2960553102122976),
    ]:
        b = r - q
        analytic = theta(is_call, S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.theta(is_call, S, K, T, r, b, v, method='central')
        assert is_close_to(numeric, analytic, 1e-4)


def test_vega():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.895452325799033),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 13.895452325799033),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.76999042895533),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 26.76999042895533),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.7352764530692),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 17.7352764530692),
    ]:
        b = r - q
        analytic = vega(S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, 1e-3)


def test_rho():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 41.58593235646499),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -5.975538868570712),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 25.08784228089154),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -22.47362894414416),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 8.182419352133445),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -44.13519899540583),
    ]:
        b = r - q
        analytic = rho(is_call, S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, 1e-4)


def test_elasticity():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 8.51357496716671),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -22.635066226567563),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 13.966967482182572),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -14.427244161400045),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 21.763120448838848),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -8.446279972615066),
    ]:
        b = r - q
        analytic = elasticity(is_call, S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.elasticity(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, 1e-4)


def test_dgamma_dvol():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0006138389948689551),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0006138389948689551),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.00338939132019472),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.00338939132019472),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.00015979636723250932),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.00015979636723250932),
    ]:
        b = r - q
        analytic = dgamma_dvol(S, K, T, r, b, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.dgamma_dvol(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, 1e-4)


def test_gammap():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.020211567019344047),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.020211567019344047),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042831984686328525),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042831984686328525),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376442324910798),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376442324910798),
    ]:
        b = r - q
        analytic = gammap(S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.gammap(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, analytic, 1e-5)


def test_vanna():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.01639625858611978),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.01639625858611978),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.002088059253458516),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.002088059253458516),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.02025315899821502),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.02025315899821502),
    ]:
        b = r - q
        analytic = vanna(S, K, T, r, b, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vanna(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, analytic, 1e-4)


def test_charm():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.23306930617480232),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.15620615104261645),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.01632708081503694),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.0931902359472228),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.2961949663176757),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.37305812144986156),
    ]:
        b = r - q
        analytic = charm(is_call, S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.charm(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, analytic, 1e-5)


def test_vegap():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 17.36931540724879),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 17.36931540724879),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 33.462488036194166),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 33.462488036194166),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 22.169095566336562),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 22.169095566336562),
    ]:
        b = r - q
        analytic = vegap(S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vegap(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, 1e-3)


def test_vomma():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0157585192593357),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0157585192593357),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.00023229659194725993),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.00023229659194725993),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.013189493867252218),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.013189493867252218),
    ]:
        b = r - q
        analytic = vomma(S, K, T, r, b, v) / 10000
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.dvega_dvol(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, 1e-4)
