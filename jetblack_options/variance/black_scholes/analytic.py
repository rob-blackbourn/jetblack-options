"""Plain Vanilla"""

from math import exp, log, sqrt

from ...distributions import CND

# The generalized Black and Scholes formula on variance form
def price(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v / 2) * T) / sqrt(v * T)
    d2 = d1 - sqrt(v * T)

    if is_call:
        return S * exp((b - r) * T) * CND(d1) - X * exp(-r * T) * CND(d2)
    else:
        return X * exp(-r * T) * CND(-d2) - S * exp((b - r) * T) * CND(-d1)
