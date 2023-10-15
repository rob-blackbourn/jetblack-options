"""Plain Vanilla"""

from math import exp, log, pi, sqrt
from typing import Literal, Callable

from ...distributions import CND, ND, CNDEV, CHIINV


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    """The generalised Black-Scholes pricing model for European options.

    The cost of carry rate (b) is:

    * b == r: for non dividend paying stocks
    * b == r - q: For dividend paying stocks where the dividend yield is q
    * b == 0: for futures options
    * b = r - rj: for currency options.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.
        cnd (Callable[[float], float], optional): The cumulative normal density
            function. Defaults to CND.

    Returns:
        float: The price of the options.
    """

    d1 = (log(S / K) + T * (b + v ** 2 / 2)) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return S * exp((b - r) * T) * cnd(d1) - K * exp(-r * T) * cnd(d2)
    else:
        return K * exp(-r * T) * cnd(-d2) - S * exp((b - r) * T) * cnd(-d1)


def delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    """The sensitivity of the open to a change in the asset price.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.
        cnd (Callable[[float], float], optional): The cumulative normal density
            function. Defaults to CND.

    Returns:
        float: the delta.
    """
    d1 = (log(S / K) + T * (b + v ** 2 / 2)) / (v * sqrt(T))

    if is_call:
        return exp((b - r) * T) * cnd(d1)
    else:
        return -exp((b - r) * T) * cnd(-d1)


def gamma(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    d1 = (log(S / K) + T * (b + v ** 2 / 2)) / (v * sqrt(T))

    return exp((b - r) * T) * nd(d1) / (S * v * sqrt(T))


def forward_delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:

    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))

    if is_call:
        return exp(-r * T) * cnd(d1)
    else:
        return exp(-r * T) * (cnd(d1) - 1)


def ddelta_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DDeltaDvol also known as vanna
    d1 = (log(S / K) + (b + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return -exp((b - r) * T) * d2 / v * nd(d1)

# DDeltaDvolDvol also known as DVannaDvol


def ddelta_dvol_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    d1 = (log(S / K) + (b + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return ddelta_dvol(S, K, T, r, b, v, nd=nd) * 1 / v * (d1 * d2 - d1 / d2 - 1)

# DdeltaDtime/Charm for the generalized Black and Scholes formula


def ddelta_dtime(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND,
        nd: Callable[[float], float] = ND
) -> float:

    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return -exp((b - r) * T) * (
            nd(d1) * (b / (v * sqrt(T)) - d2 / (2 * T)) + (b - r) * cnd(d1)
        )
    else:
        return -exp((b - r) * T) * (
            nd(d1) * (b / (v * sqrt(T)) - d2 / (2 * T)) - (b - r) * cnd(-d1)
        )


def saddle_gamma(
        K: float,
        r: float,
        b: float,
        v: float
) -> float:
    # SaddleGamma for the generalized Black and Scholes formula
    return sqrt(exp(1) / pi) * sqrt((2 * b - r) / v ** 2 + 1) / K


def dgamma_dspot(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DgammaDspot/Speed for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -gamma(S, K, T, r, b, v, nd=nd) * (1 + d1 / (v * sqrt(T))) / S


def dgamma_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DgammaDvol/Zomma for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return gamma(S, K, T, r, b, v, nd=nd) * ((d1 * d2 - 1) / v)


def dgamma_dtime(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # GGammaDtime for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return gamma(S, K, T, r, b, v, nd=nd) * (
        r - b + b * d1 / (v * sqrt(T)) +
        (1 - d1 * d2) / (2 * T)
    )


def gammap(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # GammaP for the generalized Black and Scholes formula
    return S * gamma(S, K, T, r, b, v, nd=nd) / 100


def dgammap_dspot(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DgammaPDspot/SpeedP for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -gamma(S, K, T, r, b, v, nd=nd) * (d1) / (100 * v * sqrt(T))


def dgammap_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DgammaPDvol for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S / 100 * gamma(S, K, T, r, b, v, nd=nd) * ((d1 * d2 - 1) / v)


def dgammap_dtime(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # GGammaPDtime for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        gammap(S, K, T, r, b, v, nd=nd) *
        (r - b + b * d1 / (v * sqrt(T)) + (1 - d1 * d2) / (2 * T))
    )


def vega(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # Vega for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * exp((b - r) * T) * nd(d1) * sqrt(T)


def dvega_dtime(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DvegaDtime for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        vega(S, K, T, r, b, v, nd=nd) *
        (r - b + b * d1 / (v * sqrt(T)) - (1 + d1 * d2) / (2 * T))
    )


def dvega_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DvegaDvol/Vomma for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return vega(S, K, T, r, b, v, nd=nd) * d1 * d2 / v


def dvomma_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DVommaDVol for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        dvega_dvol(S, K, T, r, b, v, nd=nd) *
        1 / v * (d1 * d2 - d1 / d2 - d2 / d1 - 1)
    )


def vegap(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # VegaP for the generalized Black and Scholes formula
    return v / 10 * vega(S, K, T, r, b, v, nd=nd)


def dvegap_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DvegaPDvol/VommaP for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return vegap(S, K, T, r, b, v, nd=nd) * d1 * d2 / v


def vega_leverage(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND,
        nd: Callable[[float], float] = ND
) -> float:
    # Vega for the generalized Black and Scholes formula
    return (
        vega(S, K, T, r, b, v, nd=nd) * v /
        price(is_call, S, K, T, r, b, v, cnd=cnd)
    )


def variance_vega(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # Variance-vega for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * exp((b - r) * T) * nd(d1) * sqrt(T) / (2 * v)


def variance_delta(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # Variance-delta for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S * exp((b - r) * T) * nd(d1) * (-d2) / (2 * v ** 2)


def variance_vomma(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # Variance-vomma for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S * exp((b - r) * T) * sqrt(T) / (4 * v ** 3) * nd(d1) * (d1 * d2 - 1)


def variance_ultima(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # Variance-ultima for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        S * exp((b - r) * T) * sqrt(T) /
        (8 * v ** 5) * nd(d1) *
        ((d1 * d2 - 1) *
         (d1 * d2 - 3) -
         (d1 ** 2 + d2 ** 2))
    )


def theta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND,
        nd: Callable[[float], float] = ND
) -> float:
    """The theta or time decay of the value of the option.

    This value is typically reported by dividing by 365 (for a one calendar day
    movement) or 252 (for a 1 trading day movement).

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        b (float): The cost of carry.
        v (float): The asset volatility.
        cnd (Callable[[float], float], optional): The cumulative density function. Defaults to CND.
        nd (Callable[[float], float], optional): The probability density function. Defaults to ND.

    Returns:
        float: The theta.
    """
    # Theta for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        p1 = -S * exp((b - r) * T) * nd(d1) * v / (2 * sqrt(T))
        p2 = (b - r) * S * exp((b - r) * T) * cnd(d1)
        p3 = r * K * exp(-r * T) * cnd(d2)
        return p1 - p2 - p3
    else:
        p1 = -S * exp((b - r) * T) * nd(d1) * v / (2 * sqrt(T))
        p2 = (b - r) * S * exp((b - r) * T) * cnd(-d1)
        p3 = r * K * exp(-r * T) * cnd(-d2)
        return p1 + p2 + p3


def theta_driftless(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # Driftless Theta for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -S * exp((b - r) * T) * nd(d1) * v / (2 * sqrt(T))


def rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    # Rho for the generalized Black and Scholes formula for all options except futures
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return T * K * exp(-r * T) * cnd(d2)
    else:
        return -T * K * exp(-r * T) * cnd(-d2)


def futures_rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    # Rho for the generalized Black and Scholes formula for Futures option
    return -T * price(is_call, S, K, T, r, 0, v, cnd=cnd)


def elasticity(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND,
) -> float:
    # Elasticity for the generalized Black and Scholes formula
    return (
        delta(is_call, S, K, T, r, b, v, cnd=cnd) * S /
        price(is_call, S, K, T, r, b, v, cnd=cnd)
    )


def carry(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    # Carry rho sensitivity for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return T * S * exp((b - r) * T) * cnd(d1)
    else:
        return -T * S * exp((b - r) * T) * cnd(-d1)


def phi(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    # Rho2/Phi for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return -T * S * exp((b - r) * T) * cnd(d1)
    else:
        return T * S * exp((b - r) * T) * cnd(-d1)


def dzeta_dvol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DZetaDvol for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return -nd(d2) * d1 / v
    else:
        return nd(d2) * d1 / v


def dzeta_dtime(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # DZetaDtime for the generalized Black and Scholes formula
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return nd(d2) * (b / (v * sqrt(T)) - d1 / (2 * T))
    else:
        return -nd(d2) * (b / (v * sqrt(T)) - d1 / (2 * T))


def break_even_probability(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    # Risk neutral break even probability for the generalized Black and Scholes formula
    if is_call:
        K = K + price(True, S, K, T, r, b, v, cnd=cnd) * exp(r * T)
        d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
        return cnd(d2)
    else:
        K = K - price(False, S, K, T, r, b, v, cnd=cnd) * exp(r * T)
        d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
        return cnd(-d2)


def strike_delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    # StrikeDelta for the generalized Black and Scholes formula
    d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return -exp(-r * T) * cnd(d2)
    else:
        return exp(-r * T) * cnd(-d2)


def risk_neutral_density(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND
) -> float:
    # Risk Neutral Density for the generalized Black and Scholes formula
    d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    return exp(-r * T) * nd(d2) / (K * v * sqrt(T))


def gamma_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        cndev: Callable[[float], float] = CNDEV
) -> float:
    # Gamma from delta
    return exp((b - r) * T) * ND(
        cndev(exp((r - b) * T) * abs(delta_))
    ) / (S * v * sqrt(T))


def gammap_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        cndev: Callable[[float], float] = CNDEV
) -> float:
    # GammaP from delta
    return S / 100 * gamma_from_delta(S, T, r, b, v, delta_, cndev=cndev)


def vega_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        delta_: float,
        *,
        cndev: Callable[[float], float] = CNDEV,
        nd: Callable[[float], float] = ND
) -> float:
    # Vega from delta
    return S * exp((b - r) * T) * sqrt(T) * nd(
        cndev(exp((r - b) * T) * abs(delta_))
    )


def vegap_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        cndev: Callable[[float], float] = CNDEV,
        nd: Callable[[float], float] = ND
) -> float:
    # VegaP from delta
    return v / 10 * vega_from_delta(S, T, r, b, delta_, cndev=cndev, nd=nd)


def strike_from_delta(
        is_call: bool,
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        cndev: Callable[[float], float] = CNDEV
) -> float:
    # Closed form solution to find strike given the delta
    if is_call:
        return S * exp(
            -cndev(delta_ * exp((r - b) * T)) *
            v * sqrt(T) + (b + v * v / 2) * T
        )
    else:
        return S * exp(
            cndev(-delta_ * exp((r - b) * T)) *
            v * sqrt(T) + (b + v * v / 2) * T
        )


def in_the_money_prob_from_delta(
        is_call: bool,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        cnd: Callable[[float], float] = CND,
        cndev: Callable[[float], float] = CNDEV
) -> float:
    # Closed form solution to find in-the-money risk-neutral probaility given the delta

    if is_call:
        return cnd(cndev(delta_ / exp((b - r) * T)) - v * sqrt(T))
    else:
        return cnd(cndev(-delta_ / exp((b - r) * T)) + v * sqrt(T))


def strike_from_in_the_money_prob(
        is_call: bool,
        S: float,
        v: float,
        T: float,
        b: float,
        in_the_money_prob: float,
        *,
        cndev: Callable[[float], float] = CNDEV
) -> float:
    # Closed form solution to find strike given the in-the-money risk neutral probability
    if is_call:
        return S * exp(
            -cndev(in_the_money_prob) * v * sqrt(T) + (b - v * v / 2) * T
        )
    else:
        return S * exp(
            cndev(in_the_money_prob) * v * sqrt(T) + (b - v * v / 2) * T
        )


def rnd_from_in_the_money_prob(
        K: float,
        T: float,
        r: float,
        v: float,
        in_the_money_prob: float,
        *,
        nd: Callable[[float], float] = ND,
        cndev: Callable[[float], float] = CNDEV
) -> float:
    # Risk Neutral Density from in-the-money probability
    return exp(-r * T) * nd(cndev(in_the_money_prob)) / (K * v * sqrt(T))


def delta_from_in_the_money_prob(
        is_call: bool,
        T: float,
        r: float,
        b: float,
        v: float,
        in_the_money_prob: float,
        *,
        cnd: Callable[[float], float] = CND,
        cndev: Callable[[float], float] = CNDEV
) -> float:
    # Closed form solution to find in-the-money risk-neutral probaility given the delta
    if is_call:
        return cnd(cndev(in_the_money_prob * exp((b - r) * T)) - v * sqrt(T))
    else:
        return -cnd(cndev(in_the_money_prob * exp((b - r) * T)) + v * sqrt(T))


def max_ddelta_dvol_asset(
        is_lower: bool,
        K: float,
        T: float,
        b: float,
        v: float
) -> float:
    # What asset price that gives maximum DdeltaDvol

    # is_lower == True gives lower asset level that gives max DdeltaDvol
    # is_lower == False gives upper asset level that gives max DdeltaDvol

    if is_lower:
        return K * exp(-b * T - v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)
    else:
        return K * exp(-b * T + v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)


def max_ddelta_dvol_strike(
        is_lower: bool,
        S: float,
        T: float,
        b: float,
        v: float
) -> float:
    # What strike price that gives maximum DdeltaDvol

    # UpperLowerFlag"l" gives lower strike level that gives max DdeltaDvol
    # UpperLowerFlag"l" gives upper strike level that gives max DdeltaDvol

    if is_lower:
        return S * exp(b * T - v * sqrt(T) * sqrt(4 + T * v * 2) / 2)
    else:
        return S * exp(b * T + v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)


def max_gamma_vega_at_X(S: float, b: float, T: float, v: float) -> float:
    # What strike price that gives maximum gamma and vega
    return S * exp((b + v * v / 2) * T)


def max_gamma_at_S(x: float, b: float, T: float, v: float) -> float:
    # What asset price that gives maximum gamma
    return x * exp((-b - 3 * v * v / 2) * T)


def max_vega_at_S(K: float, b: float, T: float, v: float) -> float:
    # What asset price that gives maximum vega
    return K * exp((-b + v * v / 2) * T)


def in_the_money_probability(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND,
) -> float:
    # Delta for the generalized Black and Scholes formula
    d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))

    if is_call:
        return cnd(d2)
    else:
        return cnd(-d2)


def delta_mirror_strike(S: float, T: float, b: float, v: float) -> float:
    # MirrorDeltaStrike, delta neutral straddle strike in the BSM formula
    return S * exp((b + v ** 2 / 2) * T)


def probability_mirror_strike(S: float, T: float, b: float, v: float) -> float:
    # MirrorProbabilityStrike, probability neutral straddle strike in the BSM formula
    return S * exp((b - v ** 2 / 2) * T)


def delta_mirror_call_put_strike(S: float, K: float, T: float, b: float, v: float) -> float:
    # MirrorProbabilityStrike, probability neutral straddle strike in the BSM formula
    return S ** 2 / K * exp((2 * b + v ** 2) * T)


def confidence_interval_volatility(
        Alfa: float,
        n: int,
        VolatilityEstimate: float,
        is_lower: bool,
        *,
        chiinv: Callable[[float, int], float] = CHIINV,
) -> float:
    # Volatility estimate confidence interval

    # UpperLower     ="L" gives the lower cofidence interval
    #                ="U" gives the upper cofidence interval
    # n: number of observations
    if is_lower:
        return VolatilityEstimate * sqrt((n - 1) / (chiinv(Alfa / 2, n - 1)))
    else:
        return VolatilityEstimate * sqrt((n - 1) / (chiinv(1 - Alfa / 2, n - 1)))


def profit_loss_std(
        TypeFlag: Literal['a', 'p'],
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        NHedges: int,
        *,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    # Profit/Loss STD for the generalized Black and Scholes formula

    if TypeFlag == "a":  # in dollars
        return (
            sqrt(pi / 4) *
            vega(S, K, T, r, b, v, nd=nd) *
            v / sqrt(NHedges)
        )
    elif TypeFlag == "p":  # in percent
        return (
            sqrt(pi / 4) *
            vega(S, K, T, r, b, v, nd=nd) *
            v / sqrt(NHedges) /
            price(is_call, S, K, T, r, b, v, cnd=cnd)
        )
    else:
        raise ValueError('invalid type flag')
