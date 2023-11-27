"""Jarrow-Rudd"""

from math import exp, nan, sqrt
from typing import Tuple


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
    """Jarrow-Rudd binomial option pricing tree.

    Args:
        is_european (bool): Tue for European, false for American.
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
    u = exp((b - v ** 2 / 2) * dT + v * sqrt(dT))
    d = exp((b - v ** 2 / 2) * dT - v * sqrt(dT))
    p = 0.5
    df = exp(-r * dT)

    option_value = [
        max(0, z * (S * u ** i * d ** (n - i) - K))
        for i in range(n+1)

    ]

    delta = gamma = theta = nan

    for j in range(n-1, -1, -1):
        for i in range(j+1):
            if is_european:
                option_value[i] = (
                    p * option_value[i + 1]
                    + (1 - p) * option_value[i]
                ) * df
            else:
                option_value[i] = max(
                    (z * (S * u ** i * d ** (j - i) - K)),
                    (
                        p * option_value[i + 1]
                        + (1 - p) * option_value[i]
                    ) * df
                )

        if j == 2:
            gamma = (
                (option_value[2] - option_value[1]) / (S * u ** 2 - S * u * d)
                - (option_value[1] - option_value[0]) /
                (S * u * d - S * d ** 2)
            ) / (0.5 * (S * u ** 2 - S * d ** 2))
            theta = option_value[1]

        if j == 1:
            delta = (option_value[1] - option_value[0]) / (S * u - S * d)

    theta = (theta - option_value[0]) / (2 * dT) / 365

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
    p, *_ = greeks(is_european, is_call, S, K, T, r, b, v, n)
    return p
