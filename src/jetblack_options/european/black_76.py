"""Black (1976) Options on futures/forwards"""

from math import exp, log, sqrt
from statistics import NormalDist

from ..implied_volatility import solve_ivol

norm = NormalDist()
cdf = norm.cdf
pdf = norm.pdf


def price(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    """Fair value of a futures/forward using Black 76.

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The price of the future.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.

    Returns:
        float: The option price.
    """

    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return exp(-r * T) * (F * cdf(d1) - K * cdf(d2))
    else:
        return exp(-r * T) * (K * cdf(-d2) - F * cdf(-d1))


def ivol(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        p: float,
        *,
        max_iterations: int = 20,
        epsilon = 1e-8
) -> float:
    """Calculate the volatility of a Black 76 option that is implied by the price.

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        p (float): The option price.
        max_iterations (int, Optional): The maximum number of iterations before
            a price is returned. Defaults to 20.
        epsilon (float, Optional): The largest acceptable error. Defaults to 1e-8.

    Returns:
        float: The implied volatility.
    """
    return solve_ivol(
        p,
        lambda v: price(is_call, F, K, T, r, v),
        max_iterations=max_iterations,
        epsilon=epsilon
    )


def delta(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    """The sensitivity of the option to a change in the asset price
    using Black 76.

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The current futures price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The volatility.

    Returns:
        float: The delta.
    """
    d1 = (log(F / K) + T * (v ** 2 / 2)) / (v * sqrt(T))
    if is_call:
        return exp(-r * T) * cdf(d1)
    else:
        return -exp(-r * T) * cdf(-d1)


def gamma(
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    """The second derivative to the change in asset price using Black 76.

    Args:
        F (float): The current futures price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The volatility.

    Returns:
        float: The gamma.
    """
    d1 = (log(F / K) + T * (v ** 2 / 2)) / (v * sqrt(T))

    return exp(-r * T) * pdf(d1) / (F * v * sqrt(T))


def theta(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    """The change in the value of the option with respect to time to expiry
    using Black 76.

    This value is typically reported by dividing by 365 (for a one calendar day
    movement) or 252 (for a 1 trading day movement).

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The current futures price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The volatility.

    Returns:
        float: The theta.
    """
    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return (
            -F * exp(-r * T) * pdf(d1) * v / (2 * sqrt(T))
            + r * F * exp(-r * T) * cdf(d1)
            - r * K * exp(-r * T) * cdf(d2)
        )
    else:
        return (
            -F * exp(-r * T) * pdf(d1) * v / (2 * sqrt(T))
            - r * F * exp(-r * T) * cdf(-d1)
            + r * K * exp(-r * T) * cdf(-d2)
        )


def vega(
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    """The sensitivity of the options price or a change in the asset volatility
    using Black 76.

    This value is typically reported by dividing by 100 (for a 1% change in
    volatility)

    Args:
        F (float): The current futures price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The volatility.

    Returns:
        float: The vega.
    """
    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    return F * exp(-r * T) * pdf(d1) * sqrt(T)


def rho(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    """The sensitivity of the option price to a change in the risk free rate
    using Black 76.

    This value is typically reported by dividing by 100 (for a 1% change
    in the risk free rate)

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The price of the future.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.

    Returns:
        float: The rho.
    """
    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return -T * exp(-r * T) * (F * cdf(d1) - K * cdf(d2))
    else:
        return -T * exp(-r * T) * (K * cdf(-d2) - F * cdf(-d1))


def vanna(
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    d1 = (log(F / K) + T * (v ** 2 / 2)) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return -exp(-r * T) * pdf(d1) * d2 / v 


def vomma(
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    d1 = (log(F / K) + T * (v ** 2 / 2)) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return F * exp(-r * T) * pdf(d1) * sqrt(T) * d1 * d2 / v
