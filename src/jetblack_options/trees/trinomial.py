"""Option valuations using a trinomial tree.
"""

from math import exp, nan, sqrt
from typing import Tuple

from ..implied_volatility import solve_ivol
from ..numeric_greeks.with_carry import NumericGreeks


def greeks(
        is_european: bool,
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        n: int
) -> Tuple[float, float, float, float]:
    """A trinomial tree options pricer returning the price and some greeks.

    Args:
        is_european (bool): True for European, false for American.
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.
        n (int): The number of the steps in the tree.

    Returns:
        Tuple[float, float, float, float]: The price, delta, gamma, theta.
    """

    z = 1 if is_call else -1

    dT = T / n
    u = exp(v * sqrt(2 * dT))
    d = exp(-v * sqrt(2 * dT))
    pu = (
        (
            exp(b * dT / 2)
            - exp(-v * sqrt(dT / 2))
        ) / (
            exp(v * sqrt(dT / 2))
            - exp(-v * sqrt(dT / 2))
        )
    ) ** 2
    pd = (
        (
            exp(v * sqrt(dT / 2))
            - exp(b * dT / 2)
        ) / (
            exp(v * sqrt(dT / 2))
            - exp(-v * sqrt(dT / 2))
        )
    ) ** 2
    pm = 1 - pu - pd
    Df = exp(-r * dT)

    option_value = [
        max(0, z * (S * u ** max(i - n, 0) * d ** max(n - i, 0) - K))
        for i in range(1 + 2*n)
    ]

    delta = gamma = theta = nan

    for j in range(n-1, -1, -1):
        for i in range(1 + j*2):

            option_value[i] = (
                pu * option_value[i + 2]
                + pm * option_value[i + 1]
                + pd * option_value[i]
            ) * Df

            if is_european:
                option_value[i] = max(
                    z * (S * u ** max(i - j, 0) * d ** max(j - i, 0) - K),
                    option_value[i]
                )

        if j == 1:
            delta = (option_value[2] - option_value[0]) / (S * u - S * d)
            gamma = (
                (option_value[2] - option_value[1]) / (S * u - S)
                - (option_value[1] - option_value[0]) / (S - S * d)
            ) / (0.5 * (S * u - S * d))
            theta = option_value[1]

    theta = (theta - option_value[0]) / dT / 365

    return option_value[0], delta, gamma, theta


def price(
        is_european: bool,
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        n: int
) -> float:
    """Calculate the price of an option using a trinomial tree.

    Args:
        is_european (bool): True for European, false for American.
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.
        n (int): The number of the steps in the tree.

    Returns:
        Tuple[float, float, float, float]: The price, delta, gamma, theta.
    """
    p, *_ = greeks(is_european, is_call, S, K, T, r, b, v, n)
    return p


def ivol(
        is_european: bool,
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        p: float,
        n: int,
        *,
        max_iterations: int = 20,
        epsilon=1e-8
) -> float:
    """Calculate the volatility of an option that is implied by the price.

    Args:
        is_european (bool): True for European, false for American.
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to expiry of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        p (float): The option price.
        n (int): The number of the steps in the tree.
        max_iterations (int, Optional): The maximum number of iterations before
            a price is returned. Defaults to 20.
        epsilon (float, Optional): The largest acceptable error. Defaults to 1e-8.

    Returns:
        float: The implied volatility.
    """
    return solve_ivol(
        p,
        lambda v: price(is_european, is_call, S, K, T, r, b, v, n),
        max_iterations=max_iterations,
        epsilon=epsilon
    )


def make_numeric_greeks(
        is_european: bool,
        is_call: bool,
        n: int
) -> NumericGreeks:
    """Make a class to generate greeks numerically using finite difference
    methods.

    Args:
        is_european (bool): True for European, false for American.
        is_call (bool): True for a call, false for a put.
        n (int): The number of the steps in the tree.

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
            b: float,
            v: float
    ) -> float:
        return price(is_european, is_call, S, K, T, r, b, v, n)

    return NumericGreeks(evaluate)
