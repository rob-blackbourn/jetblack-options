"""Leisen Reimer"""

from math import exp, log, nan, sqrt
from typing import Literal, Tuple, Union


def _sign(n: Union[float, int]) -> Literal[-1, 0, 1]:
    if n > 0:
        return 1
    elif n < 0:
        return -1
    else:
        return 0


def _odd(n: int) -> int:
    s = _sign(n)
    n = abs(n)
    if n % 2 == 0:
        n += 1
    return n * s


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
    """Leisen-Reimer binomial tree.

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

    n = _odd(n)
    z = 1 if is_call else -1

    d1 = (log(S / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    # Using Preizer-Pratt inversion method 2
    hd1 = 0.5 + _sign(d1) * (
        0.25
        - 0.25 * exp(-(d1 / (n + 1 / 3 + 0.1 / (n + 1))) ** 2 * (n + 1 / 6))
    ) ** 0.5
    hd2 = 0.5 + _sign(d2) * (
        0.25
        - 0.25 * exp(-(d2 / (n + 1 / 3 + 0.1 / (n + 1))) ** 2 * (n + 1 / 6))
    ) ** 0.5

    dT = T / n
    p = hd2
    u = exp(b * dT) * hd1 / hd2
    d = (exp(b * dT) - p * u) / (1 - p)
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
