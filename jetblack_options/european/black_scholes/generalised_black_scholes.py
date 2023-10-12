"""Generalised Black-Scholes"""

from math import exp, log, sqrt

from ...distributions import CND

class GBlackScholes:

    def __init__(self, is_call: bool, X: float, T: float) -> None:
        self.is_call = is_call
        self.X = X
        self.T = T

    
    def value(
            self,
            S: float,
            r: float,
            b: float,
            v: float
    )-> float:

        d1 = (log(S / self.X) + (b + v ** 2 / 2) * self.T) / (v * sqrt(self.T))
        d2 = d1 - v * sqrt(self.T)

        if self.is_call:
            return S * exp((b - r) * self.T) * CND(d1) - self.X * exp(-r * self.T) * CND(d2)
        else:
            return self.X * exp(-r * self.T) * CND(-d2) - S * exp((b - r) * self.T) * CND(-d1)

    def bumped_delta(self, S: float, r: float, b: float, v: float, *, dS: float = 0.01) -> float:
        return (
            self.value(S + dS, r, b, v) -
            self.value(S - dS, r, b, v)
        ) / (2 * dS)

    def delta(self, S: float, r: float, b: float, v: float) -> float:
        
        d1 = (log(S / self.X) + (b + v ** 2 / 2) * self.T) / (v * sqrt(self.T))
        if self.is_call:
            return exp((b - r) * self.T) * CND(d1)
        else:
            return -exp((b - r) * self.T) * CND(-d1)
