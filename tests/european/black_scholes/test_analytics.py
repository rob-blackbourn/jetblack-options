"""Tests for Black-Scholes European analytic options"""

from jetblack_options.european.black_scholes_merton import (
    price,
    delta,
    gamma,
    theta,
    vega,
    rho,
    elasticity,
    dgamma_dvol,
    gammap,
    ddelta_dvol,
    vegap,
    dvega_dvol,
)
from jetblack_options.numeric_greeks import NumericGreeks

from ...utils import is_close_to

def test_price():

    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 11.069546131685598, 1e-12),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.505650275001452, 1e-12),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 3.8695002999527546, 1e-12),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 2.913498834791845, 1e-12),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.7881685580252977, 1e-12),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 9.344461337871536, 1e-12),
    ]:
        b = r - q
        actual = price(is_call, S, K, T, r, b, v)
        assert is_close_to(actual, expected, 1e-12)

def test_delta():
    ng = NumericGreeks(price)

    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.8567400985874144, 1e-5),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.10404934056490878, 1e-5),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.5404518486173583, 1e-5),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.42033759053496483, 1e-5),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.17153007262292186, 1e-5),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.7892593665294013, 1e-5),
    ]:
        b = r - q
        analytic = delta(is_call, S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.delta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, precision)

def test_gamma():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018374151835767315, 1e-5),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.018374151835767315, 1e-5),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042831984686328525, 1e-5),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042831984686328525, 1e-5),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376442324910798, 1e-5),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376442324910798, 1e-5),
    ]:
        b = r - q
        analytic = gamma(S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.gamma(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, analytic, 1e-5)


def test_theta():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.006889877108078438, 1e-4),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.003993035517758725, 1e-4),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.011069047865366735, 1e-4),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.006066366408411787, 1e-4),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.006797679030347862, 1e-4),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.0008111104389378016, 1e-4),
    ]:
        b = r - q
        analytic = theta(is_call, S, K, T, r, b, v) / 365
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.theta(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, precision)


def test_vega():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13895452325799032, 1e-3),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.13895452325799032, 1e-3),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2676999042895533, 1e-3),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.2676999042895533, 1e-3),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.1773527645306925, 1e-3),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.1773527645306925, 1e-3),
    ]:
        b = r - q
        analytic = vega(S, K, T, r, b, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vega(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, precision)


def test_rho():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.41585932356464994, 1e-4),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.05975538868570709, 1e-4),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.25087842280891537, 1e-4),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.22473628944144164, 1e-4),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.08182419352133444, 1e-4),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.4413519899540583, 1e-4),
    ]:
        b = r - q
        analytic = rho(is_call, S, K, T, r, b, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.rho(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, precision)


def test_elasticity():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 8.51357496716671, 1e-4),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -22.635066226567563, 1e-4),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 13.966967482182572, 1e-4),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -14.427244161400045, 1e-4),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 21.763120448838848, 1e-4),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -8.446279972615066, 1e-4),
    ]:
        b = r - q
        analytic = elasticity(is_call, S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.elasticity(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, precision)


def test_dgamma_dvol():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0006138389948689551, 1e-4),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0006138389948689551, 1e-4),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.00338939132019472, 1e-4),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.00338939132019472, 1e-4),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.00015979636723250932, 1e-4),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, -0.00015979636723250932, 1e-4),
    ]:
        b = r - q
        analytic = dgamma_dvol(S, K, T, r, b, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.dgamma_dvol(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, precision)

def test_gammap():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.020211567019344047, 1e-5),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.020211567019344047, 1e-5),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042831984686328525, 1e-5),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.042831984686328525, 1e-5),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376442324910798, 1e-5),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.028376442324910798, 1e-5),
    ]:
        b = r - q
        analytic = gammap(S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.gammap(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, analytic, 1e-5)

def test_ddelta_dvol():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.01639625858611978, 1e-4),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, -0.01639625858611978, 1e-4),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.002088059253458516, 1e-4),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, -0.002088059253458516, 1e-4),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.02025315899821502, 1e-4),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.02025315899821502, 1e-4),
    ]:
        b = r - q
        analytic = ddelta_dvol(S, K, T, r, b, v) / 100
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.ddelta_dvol(is_call, S, K, T, r, b, v, dS=0.01)
        assert is_close_to(numeric, analytic, precision)


def test_vegap():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.17369315407248792, 1e-3),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.17369315407248792, 1e-3),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.33462488036194166, 1e-3),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.33462488036194166, 1e-3),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.22169095566336564, 1e-3),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.22169095566336564, 1e-3),
    ]:
        b = r - q
        analytic = vegap(S, K, T, r, b, v)
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.vegap(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, precision)


def test_dvega_dvol():
    ng = NumericGreeks(price)
    for is_call, S, K, r, q, T, v, expected, precision in [
        (True, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0157585192593357, 1e-4),
        (False, 110, 100, 0.1, 0.08, 6/12, 0.125, 0.0157585192593357, 1e-4),
        (True, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.00023229659194725993, 1e-4),
        (False, 100, 100, 0.1, 0.08, 6/12, 0.125, 0.00023229659194725993, 1e-4),
        (True, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.013189493867252218, 1e-4),
        (False, 100, 110, 0.1, 0.08, 6/12, 0.125, 0.013189493867252218, 1e-4),
    ]:
        b = r - q
        analytic = dvega_dvol(S, K, T, r, b, v) / 10000
        assert is_close_to(analytic, expected, 1e-12)

        numeric = ng.dvega_dvol(is_call, S, K, T, r, b, v)
        assert is_close_to(numeric, analytic, precision)
