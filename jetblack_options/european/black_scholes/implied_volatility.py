"""Black Scholes Merton"""

from typing import Callable

from ...distributions import CND

from .analytic import price

def ivol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        p: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:

    vLow = 0.005
    vHigh = 4
    epsilon = 0.00000001
    cLow = price(is_call, S, K, T, r, b, vLow, cnd=cnd)
    cHigh = price(is_call, S, K, T, r, b, vHigh, cnd=cnd)
    N = 0
    vi = vLow + (p - cLow) * (vHigh - vLow) / (cHigh - cLow)
    while abs(p - price(is_call, S, K, T, r, b, vi, cnd=cnd)) > epsilon:
        N = N + 1
        if N > 20:
            break
        
        if price(is_call, S, K, T, r, b, vi, cnd=cnd) < p:
            vLow = vi
        else:
            vHigh = vi

        cLow = price(is_call, S, K, T, r, b, vLow, cnd=cnd)
        cHigh = price(is_call, S, K, T, r, b, vHigh, cnd=cnd)
        vi = vLow + (p - cLow) * (vHigh - vLow) / (cHigh - cLow)

    return vi
