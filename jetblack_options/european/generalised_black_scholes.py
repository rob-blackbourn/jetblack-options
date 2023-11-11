"""Black-Scholes-Merton options pricing formulae.

This is the "generalised" version using "cost of carry" (variable b).

The cost of carry rate (b) is:

* b == r: for non dividend paying stocks
* b == r - q: For dividend paying stocks where the dividend yield is q
* b == 0: for futures options
* b = r - rj: for currency options.
"""

from math import exp, log, pi, sqrt
from typing import Literal, Callable

from ..distributions import CDF, PDF, INV_CDF, CHIINV

from ..implied_volatility import solve_ivol


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    """The fair value of a European option, using Black-Scholes-Merton.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.
        cdf (Callable[[float], float], optional): The cumulative probability
            distribution function. Defaults to CDF.

    Returns:
        float: The price of the options.
    """

    d1 = (log(S / K) + T * (b + v ** 2 / 2)) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return S * exp((b - r) * T) * cdf(d1) - K * exp(-r * T) * cdf(d2)
    else:
        return K * exp(-r * T) * cdf(-d2) - S * exp((b - r) * T) * cdf(-d1)


def ivol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        p: float,
        *,
        cdf: Callable[[float], float] = CDF,
        max_iterations: int = 20,
        epsilon = 1e-8
) -> float:
    """Calculate the volatility of an option that is implied by the price.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        p (float): The option price.
        cdf (Callable[[float], float], optional): The cumulative probability
            distribution function. Defaults to CDF.
        max_iterations (int, Optional): The maximum number of iterations before
            a price is returned. Defaults to 20.
        epsilon (float, Optional): The largest acceptable error. Defaults to 1e-8.

    Returns:
        float: The implied volatility.
    """
    return solve_ivol(
        p,
        lambda v: price(is_call, S, K, T, r, b, v, cdf=cdf),
        max_iterations=max_iterations,
        epsilon=epsilon
    )


def delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
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
        cdf (Callable[[float], float], optional): The cumulative probability
            distribution function. Defaults to CDF.

    Returns:
        float: the delta.
    """
    d1 = (log(S / K) + T * (b + v ** 2 / 2)) / (v * sqrt(T))

    if is_call:
        return exp((b - r) * T) * cdf(d1)
    else:
        return -exp((b - r) * T) * cdf(-d1)


def gamma(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    """The second derivative to the change in the asset price.

    Args:
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.
        pdf (Callable[[float], float], optional): The probability distribution
            function. Defaults to PDF.

    Returns:
        float: The gamma.
    """
    d1 = (log(S / K) + T * (b + v ** 2 / 2)) / (v * sqrt(T))

    return exp((b - r) * T) * pdf(d1) / (S * v * sqrt(T))


def theta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF,
        pdf: Callable[[float], float] = PDF
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
        cdf (Callable[[float], float], optional): The cumulative density function. Defaults to CDF.
        pdf (Callable[[float], float], optional): The probability density function. Defaults to PDF.

    Returns:
        float: The theta.
    """
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        p1 = -S * exp((b - r) * T) * pdf(d1) * v / (2 * sqrt(T))
        p2 = (b - r) * S * exp((b - r) * T) * cdf(d1)
        p3 = r * K * exp(-r * T) * cdf(d2)
        return p1 - p2 - p3
    else:
        p1 = -S * exp((b - r) * T) * pdf(d1) * v / (2 * sqrt(T))
        p2 = (b - r) * S * exp((b - r) * T) * cdf(-d1)
        p3 = r * K * exp(-r * T) * cdf(-d2)
        return p1 + p2 + p3


def vega(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    """The sensitivity of the options price or a change in the asset volatility.

    Args:
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.
        pdf (Callable[[float], float], optional): The probability distribution
            function. Defaults to PDF.

    Returns:
        float: The vega
    """
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * exp((b - r) * T) * pdf(d1) * sqrt(T)


def rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    """The sensitivity of the option price to the risk free rate.

    Useful for all options except futures options which should use
    futures_rho.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        b (float): The cost of carry.
        v (float): The asset volatility.
        cdf (Callable[[float], float], optional): The cumulative density function. Defaults to CDF.

    Returns:
        float: The rho.
    """
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return T * K * exp(-r * T) * cdf(d2)
    else:
        return -T * K * exp(-r * T) * cdf(-d2)


def carry(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    """Sensitivity to the cost of carry.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        b (float): The cost of carry.
        v (float): The asset volatility.
        cdf (Callable[[float], float], optional): The cumulative density function. Defaults to CDF.
        pdf (Callable[[float], float], optional): The probability density function. Defaults to PDF.

    Returns:
        float: The carry.
    """
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return T * S * exp((b - r) * T) * cdf(d1)
    else:
        return -T * S * exp((b - r) * T) * cdf(-d1)


def elasticity(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF,
) -> float:
    """The option elasticity.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        b (float): The cost of carry.
        v (float): The asset volatility.
        cdf (Callable[[float], float], optional): The cumulative density function. Defaults to CDF.

    Returns:
        float: The elasticity.
    """
    return (
        delta(is_call, S, K, T, r, b, v, cdf=cdf) * S /
        price(is_call, S, K, T, r, b, v, cdf=cdf)
    )


def gammap(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    return S * gamma(S, K, T, r, b, v, pdf=pdf) / 100


def vegap(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    return v / 10 * vega(S, K, T, r, b, v, pdf=pdf)


def forward_delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:

    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))

    if is_call:
        return exp(-r * T) * cdf(d1)
    else:
        return exp(-r * T) * (cdf(d1) - 1)


def vanna(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    # Also known as DdeltaDvol.
    d1 = (log(S / K) + (b + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return -exp((b - r) * T) * d2 / v * pdf(d1)


def ddelta_dvol_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    # Also known as DVannaDvol
    d1 = (log(S / K) + (b + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return vanna(S, K, T, r, b, v, pdf=pdf) * 1 / v * (d1 * d2 - d1 / d2 - 1)


def charm(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF,
        pdf: Callable[[float], float] = PDF
) -> float:
    # Also known as DdeltaDtime

    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return -exp((b - r) * T) * (
            pdf(d1) * (b / (v * sqrt(T)) - d2 / (2 * T)) + (b - r) * cdf(d1)
        )
    else:
        return -exp((b - r) * T) * (
            pdf(d1) * (b / (v * sqrt(T)) - d2 / (2 * T)) - (b - r) * cdf(-d1)
        )


def saddle_gamma(
        K: float,
        r: float,
        b: float,
        v: float
) -> float:
    return sqrt(exp(1) / pi) * sqrt((2 * b - r) / v ** 2 + 1) / K


def dgamma_dspot(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    # Also known as Speed
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -gamma(S, K, T, r, b, v, pdf=pdf) * (1 + d1 / (v * sqrt(T))) / S


def dgamma_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    # Also known as zomma.
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return gamma(S, K, T, r, b, v, pdf=pdf) * ((d1 * d2 - 1) / v)


def dgamma_dtime(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return gamma(S, K, T, r, b, v, pdf=pdf) * (
        r - b + b * d1 / (v * sqrt(T)) +
        (1 - d1 * d2) / (2 * T)
    )


def dgammap_dspot(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    # Also known as SpeedP.
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -gamma(S, K, T, r, b, v, pdf=pdf) * (d1) / (100 * v * sqrt(T))


def dgammap_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S / 100 * gamma(S, K, T, r, b, v, pdf=pdf) * ((d1 * d2 - 1) / v)


def dgammap_dtime(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        gammap(S, K, T, r, b, v, pdf=pdf) *
        (r - b + b * d1 / (v * sqrt(T)) + (1 - d1 * d2) / (2 * T))
    )


def dvega_dtime(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        vega(S, K, T, r, b, v, pdf=pdf) *
        (r - b + b * d1 / (v * sqrt(T)) - (1 + d1 * d2) / (2 * T))
    )


def vomma(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    # Also known as DvegaDvol
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return vega(S, K, T, r, b, v, pdf=pdf) * d1 * d2 / v


def dvomma_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        vomma(S, K, T, r, b, v, pdf=pdf) *
        1 / v * (d1 * d2 - d1 / d2 - d2 / d1 - 1)
    )


def dvegap_dvol(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    # Also known as VommaP.
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return vegap(S, K, T, r, b, v, pdf=pdf) * d1 * d2 / v


def vega_leverage(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF,
        pdf: Callable[[float], float] = PDF
) -> float:
    return (
        vega(S, K, T, r, b, v, pdf=pdf) * v /
        price(is_call, S, K, T, r, b, v, cdf=cdf)
    )


def variance_vega(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * exp((b - r) * T) * pdf(d1) * sqrt(T) / (2 * v)


def variance_delta(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S * exp((b - r) * T) * pdf(d1) * (-d2) / (2 * v ** 2)


def variance_vomma(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S * exp((b - r) * T) * sqrt(T) / (4 * v ** 3) * pdf(d1) * (d1 * d2 - 1)


def variance_ultima(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        S * exp((b - r) * T) * sqrt(T) /
        (8 * v ** 5) * pdf(d1) *
        ((d1 * d2 - 1) *
         (d1 * d2 - 3) -
         (d1 ** 2 + d2 ** 2))
    )


def theta_driftless(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -S * exp((b - r) * T) * pdf(d1) * v / (2 * sqrt(T))


def futures_rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    return -T * price(is_call, S, K, T, r, 0, v, cdf=cdf)


def phi(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    # Also known as rho2.
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return -T * S * exp((b - r) * T) * cdf(d1)
    else:
        return T * S * exp((b - r) * T) * cdf(-d1)


def dzeta_dvol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return -pdf(d2) * d1 / v
    else:
        return pdf(d2) * d1 / v


def dzeta_dtime(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return pdf(d2) * (b / (v * sqrt(T)) - d1 / (2 * T))
    else:
        return -pdf(d2) * (b / (v * sqrt(T)) - d1 / (2 * T))


def break_even_probability(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    # Risk neutral break even probability.
    if is_call:
        K = K + price(True, S, K, T, r, b, v, cdf=cdf) * exp(r * T)
        d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
        return cdf(d2)
    else:
        K = K - price(False, S, K, T, r, b, v, cdf=cdf) * exp(r * T)
        d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
        return cdf(-d2)


def strike_delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return -exp(-r * T) * cdf(d2)
    else:
        return exp(-r * T) * cdf(-d2)


def risk_neutral_density(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF
) -> float:
    d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    return exp(-r * T) * pdf(d2) / (K * v * sqrt(T))


def gamma_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        inv_cdf: Callable[[float], float] = INV_CDF
) -> float:
    return exp((b - r) * T) * PDF(
        inv_cdf(exp((r - b) * T) * abs(delta_))
    ) / (S * v * sqrt(T))


def gammap_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        inv_cdf: Callable[[float], float] = INV_CDF
) -> float:
    return S / 100 * gamma_from_delta(S, T, r, b, v, delta_, inv_cdf=inv_cdf)


def vega_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        delta_: float,
        *,
        inv_cdf: Callable[[float], float] = INV_CDF,
        pdf: Callable[[float], float] = PDF
) -> float:
    return S * exp((b - r) * T) * sqrt(T) * pdf(
        inv_cdf(exp((r - b) * T) * abs(delta_))
    )


def vegap_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        inv_cdf: Callable[[float], float] = INV_CDF,
        pdf: Callable[[float], float] = PDF
) -> float:
    return v / 10 * vega_from_delta(S, T, r, b, delta_, inv_cdf=inv_cdf, pdf=pdf)


def strike_from_delta(
        is_call: bool,
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta_: float,
        *,
        inv_cdf: Callable[[float], float] = INV_CDF
) -> float:
    if is_call:
        return S * exp(
            -inv_cdf(delta_ * exp((r - b) * T)) *
            v * sqrt(T) + (b + v * v / 2) * T
        )
    else:
        return S * exp(
            inv_cdf(-delta_ * exp((r - b) * T)) *
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
        cdf: Callable[[float], float] = CDF,
        inv_cdf: Callable[[float], float] = INV_CDF
) -> float:
    if is_call:
        return cdf(inv_cdf(delta_ / exp((b - r) * T)) - v * sqrt(T))
    else:
        return cdf(inv_cdf(-delta_ / exp((b - r) * T)) + v * sqrt(T))


def strike_from_in_the_money_prob(
        is_call: bool,
        S: float,
        v: float,
        T: float,
        b: float,
        in_the_money_prob: float,
        *,
        inv_cdf: Callable[[float], float] = INV_CDF
) -> float:
    if is_call:
        return S * exp(
            -inv_cdf(in_the_money_prob) * v * sqrt(T) + (b - v * v / 2) * T
        )
    else:
        return S * exp(
            inv_cdf(in_the_money_prob) * v * sqrt(T) + (b - v * v / 2) * T
        )


def rnd_from_in_the_money_prob(
        K: float,
        T: float,
        r: float,
        v: float,
        in_the_money_prob: float,
        *,
        pdf: Callable[[float], float] = PDF,
        inv_cdf: Callable[[float], float] = INV_CDF
) -> float:
    return exp(-r * T) * pdf(inv_cdf(in_the_money_prob)) / (K * v * sqrt(T))


def delta_from_in_the_money_prob(
        is_call: bool,
        T: float,
        r: float,
        b: float,
        v: float,
        in_the_money_prob: float,
        *,
        cdf: Callable[[float], float] = CDF,
        inv_cdf: Callable[[float], float] = INV_CDF
) -> float:
    if is_call:
        return cdf(inv_cdf(in_the_money_prob * exp((b - r) * T)) - v * sqrt(T))
    else:
        return -cdf(inv_cdf(in_the_money_prob * exp((b - r) * T)) + v * sqrt(T))


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

    # is_lower == True gives lower strike level that gives max DdeltaDvol
    # is_lower == False gives upper strike level that gives max DdeltaDvol

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
        cdf: Callable[[float], float] = CDF,
) -> float:
    d2 = (log(S / K) + (b - v ** 2 / 2) * T) / (v * sqrt(T))

    if is_call:
        return cdf(d2)
    else:
        return cdf(-d2)


def delta_mirror_strike(S: float, T: float, b: float, v: float) -> float:
    return S * exp((b + v ** 2 / 2) * T)


def probability_mirror_strike(S: float, T: float, b: float, v: float) -> float:
    return S * exp((b - v ** 2 / 2) * T)


def delta_mirror_call_put_strike(S: float, K: float, T: float, b: float, v: float) -> float:
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

    # is_lower == True gives the lower cofidence interval
    # is_lower == False gives the upper cofidence interval
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
        pdf: Callable[[float], float] = PDF,
        cdf: Callable[[float], float] = CDF
) -> float:

    if TypeFlag == "a":  # in dollars
        return (
            sqrt(pi / 4) *
            vega(S, K, T, r, b, v, pdf=pdf) *
            v / sqrt(NHedges)
        )
    elif TypeFlag == "p":  # in percent
        return (
            sqrt(pi / 4) *
            vega(S, K, T, r, b, v, pdf=pdf) *
            v / sqrt(NHedges) /
            price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
    else:
        raise ValueError('invalid type flag')