"""Optional valuation with a European binomial implementation.
"""

from math import comb, exp, log, sqrt

from ..implied_volatility import solve_ivol
from ..numeric_greeks.with_carry import NumericGreeks


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        n: int
) -> float:
    """A European binomial option pricing tree.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.
        n (int): The number of the steps in the tree.

    Returns:
        float: The price of the option.
    """

    dt = T / n
    u = exp(v * sqrt(dt))
    d = 1 / u
    a = exp(b * dt)
    p = (a - d) / (u - d)
    A = int(log(K / (S * d ** n)) / log(u / d)) + 1

    sum = 0
    if is_call:
        for j in range(A, n+1):
            sum += comb(n, j) * p ** j * (1 - p) ** (n - j) * (
                S * u ** j * d ** (n - j) - K
            )
    else:
        for j in range(A):
            sum += comb(n, j) * p ** j * (1 - p) ** (n - j) * (
                K - S * u ** j * d ** (n - j)
            )

    return exp(-r * T) * sum


def ivol(
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
        lambda v: price(is_call, S, K, T, r, b, v, n),
        max_iterations=max_iterations,
        epsilon=epsilon
    )


def make_numeric_greeks(is_call: bool, n: int) -> NumericGreeks:
    """_summary_

    Args:
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
        return price(is_call, S, K, T, r, b, v, n)

    return NumericGreeks(evaluate)
