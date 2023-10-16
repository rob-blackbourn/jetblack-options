"""European binomial"""

from math import comb, exp, log, sqrt

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
