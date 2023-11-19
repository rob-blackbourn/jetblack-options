"""Black-Scholes 1973"""

from math import exp, log, sqrt
from statistics import NormalDist

from ..implied_volatility import solve_ivol

norm = NormalDist()
cdf = norm.cdf
pdf = norm.pdf


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    """Black-Scholes for a non-dividend paying stock.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.

    Returns:
        float: The price of the option.
    """
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return S * cdf(d1) - K * exp(-r*T) * cdf(d2)
    else:
        return K * exp(-r * T) * cdf(-d2) - S * cdf(-d1)


def ivol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        p: float,
        *,
        max_iterations: int = 35,
        epsilon=1e-8
) -> float:
    """Calculate the volatility of a Black-Scholes 73 option that is implied by
    the price.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        p (float): The option price.
        max_iterations (int, Optional): The maximum number of iterations before
            a price is returned. Defaults to 35.
        epsilon (float, Optional): The largest acceptable error. Defaults to 1e-8.

    Returns:
        float: The implied volatility.
    """
    return solve_ivol(
        p,
        lambda v: price(is_call, S, K, T, r, v),
        max_iterations=max_iterations,
        epsilon=epsilon
    )


def delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return cdf(d1)
    else:
        return -cdf(-d1)

def gamma (
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    "Calculates option gamma"
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    return pdf(d1) / (S * v * sqrt(T))

def theta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    "Calculates option theta"
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    
    if is_call:
        return - (
            (S * pdf(d1) * v) / (2 * sqrt(T))
        ) - r * K * exp(-r * T) * cdf(d2)

    else:
        return - (
            (S * pdf(d1) * v) / (2 * sqrt(T))
        ) + r * K * exp(-r * T) * cdf(-d2)


def vega (
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * sqrt(T) * pdf(d1)


def rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    
    if is_call:
        return K * T * exp(-r * T) * cdf(d2)
    else:
        return -K * T * exp(-r * T) * cdf(-d2)


def vanna(
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    # Also known as DdeltaDvol.
    d1 = (log(S / K) + (r + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return -d2 * pdf(d1) / v

def charm(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    # Also known as DdeltaDtime

    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return -pdf(d1) * (r / (v * sqrt(T)) - d2 / (2 * T))
    else:
        return -pdf(d1) * (r / (v * sqrt(T)) - d2 / (2 * T))


def vomma(
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    # Also known as DvegaDvol
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return vega(S, K, T, r, v) * d1 * d2 / v
