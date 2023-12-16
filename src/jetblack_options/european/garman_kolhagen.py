"""Garman and Kolhagen (1983) Currency options"""

from math import exp, log, sqrt
from statistics import NormalDist

from ..implied_volatility import solve_ivol
from ..numeric_greeks.with_dividend_yield import NumericGreeks

norm = NormalDist()
cdf = norm.cdf
pdf = norm.pdf
inv_cdf = norm.inv_cdf


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        rf: float,
        v: float,
) -> float:
    """Garman and Kolhagen (1983) Currency options.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate of the base currency.
        rf (float): The risk free rate of the quote currency.
        v (float): The asset volatility.

    Returns:
        float: _description_
    """
    # Garman and Kolhagen (1983) Currency options

    d1 = (log(S / K) + (r - rf + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return S * exp(-rf * T) * cdf(d1) - K * exp(-r * T) * cdf(d2)
    else:
        return K * exp(-r * T) * cdf(-d2) - S * exp(-rf * T) * cdf(-d1)


def ivol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        rf: float,
        p: float,
        *,
        max_iterations: int = 20,
        epsilon=1e-8
) -> float:
    """Calculate the volatility of an option that is implied by the price.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to expiry of the option in years.
        r (float): The risk free rate of the base currency.
        rf (float): The risk free rate of the quote currency.
        p (float): The option price.
        max_iterations (int, Optional): The maximum number of iterations before
            a price is returned. Defaults to 20.
        epsilon (float, Optional): The largest acceptable error. Defaults to 1e-8.

    Returns:
        float: The implied volatility.
    """
    return solve_ivol(
        p,
        lambda v: price(is_call, S, K, T, r, rf, v),
        max_iterations=max_iterations,
        epsilon=epsilon
    )


def make_numeric_greeks(is_call: bool) -> NumericGreeks:
    """Make a class to generate greeks numerically using finite difference
    methods.

    Args:
        is_call (bool): If true the options is a call;  otherwise it is a put.

    Returns:
        NumericGreeks: A class which can generate Greeks using finite difference
            methods.
    """
    # Normalize the price function to match that required by the finite
    # difference methods.
    def evaluate(
            S: float,
            K: float,
            T: float,
            r: float,
            rf: float,
            v: float
    ) -> float:
        return price(is_call, S, K, T, r, rf, v)

    return NumericGreeks(evaluate)
