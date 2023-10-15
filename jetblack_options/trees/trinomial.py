from math import exp, log, nan, sqrt
from typing import Tuple

def trinomial(
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

    # Dim OptionValue() As Double
    # ReDim OptionValue(0 To n * 2 + 1)

    # Dim ReturnValue(3) As Double
    # Dim dt As Double, u As Double, d As Double
    # Dim pu As Double, pd As Double, pm As Double, Df As Double
    # Dim i As Long, j As Long, z As Integer

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
        for i in range(1+ j*2):
        
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
