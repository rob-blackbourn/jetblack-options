"""Cox-Ross-Rubinstein"""

from math import exp, nan, sqrt
from typing import Tuple

# Cox-Ross-Rubinstein binomial tree
def crr_binomial(
        is_european: bool,
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        n: int
) -> Tuple[float, float, float, float]:

    z = 1 if is_call else -1
    
    dT = T / n
    u = exp(v * sqrt(dT))
    d = 1 / u
    a = exp(b * dT)
    p = (a - d) / (u - d)
    df = exp(-r * dT)
    
    option_value = [
        max(0, z * (S * u ** i * d ** (n - i) - X))
        for i in range(n+1)
    ]

    delta = gamma = theta = nan

    for j in range(n-1, -1, -1):
        for i in range(j+1):
            if is_european:
                option_value[i] = (
                    p * option_value[i + 1] +
                    (1 - p) * option_value[i]
                ) * df
            else:
                option_value[i] = max(
                    (z * (S * u ** i * d ** (j - i) - X)),
                    (
                        p * option_value[i + 1]
                        + (1 - p) * option_value[i]
                    ) * df
                )

        if j == 2:
            gamma = (
                (option_value[2] - option_value[1]) / (S * u ** 2 - S)
                - (option_value[1] - option_value[0]) / (S - S * d ** 2)
            ) / (0.5 * (S * u ** 2 - S * d ** 2))
            theta = option_value[1]
            
        if j == 1:
            delta = (option_value[1] - option_value[0]) / (S * u - S * d)

    theta = (theta - option_value[0]) / (2 * dT) / 365

    return option_value[0], delta, gamma, theta
