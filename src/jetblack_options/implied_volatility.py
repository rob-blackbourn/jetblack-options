"""implied volatility"""

from typing import Callable


def solve_ivol(
        p: float,
        price: Callable[[float], float],
        *,
        max_iterations: int = 20,
        epsilon=1e-8
) -> float:

    v_lo = 0.005
    v_hi = 4.0
    p_lo = price(v_lo)
    p_hi = price(v_hi)

    n = 0
    v = v_lo + (p - p_lo) * (v_hi - v_lo) / (p_hi - p_lo)
    p1 = price(v)
    while abs(p - p1) > epsilon and n < max_iterations:
        n += 1

        if p1 < p:
            v_lo = v
        else:
            v_hi = v

        p_lo = price(v_lo)
        p_hi = price(v_hi)
        v = v_lo + (p - p_lo) * (v_hi - v_lo) / (p_hi - p_lo)
        p1 = price(v)

    return v
