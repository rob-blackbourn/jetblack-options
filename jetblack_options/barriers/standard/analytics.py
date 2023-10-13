"""Standard barrier option - analytic"""

from math import exp, log, sqrt

from ...distributions import CND


def price(
        is_call: bool,
        is_up: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        k: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    # Standard barrier options

    mu = (b - v ** 2 / 2) / v ** 2
    lambda_ = sqrt(mu ** 2 + 2 * r / v ** 2)
    X1 = log(S / X) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T)
    X2 = log(S / H) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T)
    y1 = log(H ** 2 / (S * X)) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T)
    y2 = log(H / S) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T)
    z = log(H / S) / (v * sqrt(T)) + lambda_ * v * sqrt(T)
    
    if is_call:
        if is_up:
            eta = -1
            phi = 1
        else:
            eta = 1
            phi = 1
    else:
        if is_up:
            eta = -1
            phi = -1
        else:
            eta = 1
            phi = -1
    
    f1 = phi * S * exp((b - r) * T) * CND(phi * X1) - phi * X * exp(-r * T) * CND(phi * X1 - phi * v * sqrt(T))
    f2 = phi * S * exp((b - r) * T) * CND(phi * X2) - phi * X * exp(-r * T) * CND(phi * X2 - phi * v * sqrt(T))
    f3 = phi * S * exp((b - r) * T) * (H / S) ** (2 * (mu + 1)) * CND(eta * y1) - phi * X * exp(-r * T) * (H / S) ** (2 * mu) * CND(eta * y1 - eta * v * sqrt(T))
    f4 = phi * S * exp((b - r) * T) * (H / S) ** (2 * (mu + 1)) * CND(eta * y2) - phi * X * exp(-r * T) * (H / S) ** (2 * mu) * CND(eta * y2 - eta * v * sqrt(T))
    f5 = k * exp(-r * T) * (CND(eta * X2 - eta * v * sqrt(T)) - (H / S) ** (2 * mu) * CND(eta * y2 - eta * v * sqrt(T)))
    f6 = k * ((H / S) ** (mu + lambda_) * CND(eta * z) + (H / S) ** (mu - lambda_) * CND(eta * z - 2 * eta * lambda_ * v * sqrt(T)))
    
    if X > H:
        if is_call:
            if is_up:
                if is_in:
                    return f1 + f5
                else:
                    return f6
            else:
                if is_in:
                    return f3 + f5
                else:
                    return f1 - f3 + f6
        else:
            if is_up:
                if is_in:
                    return f1 - f2 + f4 + f5
                else:
                    return f2 - f4 + f6
            else:
                if is_in:
                    return f2 - f3 + f4 + f5
                else:
                    return f1 - f2 + f3 - f4 + f6
    elif X < H:
        if is_call:
            if is_up:
                if is_in:
                    return f2 - f3 + f4 + f5
                else:
                    return f1 - f2 + f3 - f4 + f6
            else:
                if is_in:
                    return f1 - f2 + f4 + f5
                else:
                    return f2 + f6 - f4
        else:
            if is_up:
                if is_in:
                    return f3 + f5
                else:
                    return f1 - f3 + f6
            else:
                if is_in:
                    return f1 + f5
                else:
                    return f6
    else:
        raise ValueError("no solution")
