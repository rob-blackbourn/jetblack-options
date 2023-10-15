"""Asian options"""

from math import exp, log, sqrt
from typing import Callable

from ...distributions import CND
from ...european.black_scholes.analytic import price as bs_price


def price(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CND
) -> float:

    z = 1 if is_call else -1

    dt = (T - t1) / (n - 1)
    
    if b == 0:
        EA = S
    else:
        EA = S / n * exp(b * t1) * (1 - exp(b * dt * n)) / (1 - exp(b * dt))

    if m > 0:
        if SA > n / m * K:
            # Exercise is certain for call, put must be out-of-the-money:
            if not is_call:
                return 0
            else:
                SA = SA * m / n + EA * (n - m) / n
                return (SA - K) * exp(-r * T)

    if m == n - 1:
        # Only one fix left use Black-Scholes weighted with time
        K = n * K - (n - 1) * SA
        return bs_price(is_call, S, K, T, r, b, v, cdf=cdf) * 1 / n

    if m > 0:
        K = n / (n - m) * K - m / (n - m) * SA

    vx = v * sqrt(t1 + dt * (n - 1) * (2 * n - 1) / (6 * n))
    my = log(S) + (b - v * v * 0.5) * (t1 + (n - 1) * dt / 2)

    sum1 = 0
    for i in range(1, 1+n):
    
        ti = dt * i + t1 - dt
        vi = v * sqrt(t1 + (i - 1) * dt)
        vxi = v * v * (t1 + dt * ((i - 1) - i * (i - 1) / (2 * n)))
        myi = log(S) + (b - v * v * 0.5) * ti
        sum1 = sum1 + exp(myi + vxi / (vx * vx) * (log(K) - my) + (vi * vi - vxi * vxi / (vx * vx)) * 0.5)
    Km = 2 * K - 1 / n * sum1
    sum2 = 0

    for i in range(1, 1+n):
    
        ti = dt * i + t1 - dt
        vi = v * sqrt(t1 + (i - 1) * dt)
        vxi = v * v * (t1 + dt * ((i - 1) - i * (i - 1) / (2 * n)))
        myi = log(S) + (b - v * v * 0.5) * ti
        sum2 = sum2 + exp(myi + vi * vi * 0.5) * cdf(z * ((my - log(Km)) / vx + vxi / vx))

    return exp(-r * T) * z * (1 / n * sum2 - K * cdf(z * (my - log(Km)) / vx)) * (n - m) / n
