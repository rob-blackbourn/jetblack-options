"""Perpetual analytic solutions"""

from math import sqrt

def price(
        is_call: bool,
        S: float,
        K: float,
        r: float,
        b: float,
        v: float
) -> float:

    y1 = 1 / 2 - b / v ** 2 + sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
    y2 = 1 / 2 - b / v ** 2 - sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
    if is_call:
        return K / (y1 - 1) * ((y1 - 1) / y1 * S / K) ** y1
    else:
        return K / (1 - y2) * ((y2 - 1) / y2 * S / K) ** y2
