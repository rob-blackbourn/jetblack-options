"""Plain Vanilla"""

from math import exp, log, sqrt
from typing import Literal

from ..distributions import CND


# Black (1976) Options on futures/forwards
def Black76(CallPutFlag: Literal['c', 'p'], F: float, X: float, T: float, r: float, v: float) -> float:

    d1 = (log(F / X) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if CallPutFlag == "c":
        return exp(-r * T) * (F * CND(d1) - X * CND(d2))
    elif CallPutFlag == "p":
        return exp(-r * T) * (X * CND(-d2) - F * CND(-d1))
    else:
        raise ValueError("invalid call put flag")
