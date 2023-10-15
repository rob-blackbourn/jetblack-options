"""Black-Scholes variance numeric solutions"""

from .analytic import price

def delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
            price(is_call, S + dS, K, T, r, b, v) -
            price(is_call, S - dS, K, T, r, b, v)
    ) / (2 * dS)

def gamma(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, K, T, r, b, v) -
        2 * price(is_call, S, K, T, r, b, v) +
        price(is_call, S - dS, K, T, r, b, v)
    ) / dS ** 2

def theta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365.0
) -> float:
    if T <= dT:
        return (
            price(is_call, S, K, 0.00001, r, b, v) -
            price(is_call, S, K, T, r, b, v)
        )
    else:
        return (
            price(is_call, S, K, T - dT, r, b, v) -
            price(is_call, S, K, T, r, b, v)
        )

def vega(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01
) -> float:
    return (
        price(is_call, S, K, T, r, b, v + dv) -
        price(is_call, S, K, T, r, b, v - dv)
    ) / 2


def rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01
) -> float:
    return (
        price(is_call, S, K, T, r + dr, b + dr, v) -
        price(is_call, S, K, T, r - dr, b - dr, v)
    ) / 2

def carry(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01
) -> float:
    return (
        price(is_call, S, K, T, r, b + db, v) -
        price(is_call, S, K, T, r, b - db, v)
    ) / 2


def elasticity(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, K, T, r, b, v) -
        price(is_call, S - dS, K, T, r, b, v)
    ) / (2 * dS) * S / price(is_call, S, K, T, r, b, v)

def speed(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + 2 * dS, K, T, r, b, v) -
        3 * price(is_call, S + dS, K, T, r, b, v) +
        3 * price(is_call, S, K, T, r, b, v) -
        price(is_call, S - dS, K, T, r, b, v)
    ) / dS ** 3

def gammap(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return S / 100 * (
        price(is_call, S + dS, K, T, r, b, v) -
        2 * price(is_call, S, K, T, r, b, v) +
        price(is_call, S - dS, K, T, r, b, v)
    ) / dS ** 2


def vegap(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01
) -> float:
    return v / 0.1 * (
        price(is_call, S, K, T, r, b, v + dv) -
        price(is_call, S, K, T, r, b, v - dv)
    ) / 2

def ddelta_dvariance(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01
) -> float:
    return 1 / (4 * dS * 0.01) * (
        price(is_call, S + dS, K, T, r, b, v + dv) -
        price(is_call, S + dS, K, T, r, b, v - dv) -
        price(is_call, S - dS, K, T, r, b, v + dv) +
        price(is_call, S - dS, K, T, r, b, v - dv)
    ) / 100

def dgamma_dvariance(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, K, T, r, b, v + dv) -
        2 * price(is_call, S, K, T, r, b, v + dv) +
        price(is_call, S - dS, K, T, r, b, v + dv) -
        price(is_call, S + dS, K, T, r, b, v - dv) +
        2 * price(is_call, S, K, T, r, b, v - dv) -
        price(is_call, S - dS, K, T, r, b, v - dv)
    ) / (2 * dv * dS ** 2) / 100

def dvega_variance(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01
) -> float:
    return (
        price(is_call, S, K, T, r, b, v + dv) -
        2 * price(is_call, S, K, T, r, b, v) +
        price(is_call, S, K, T, r, b, v - dv)
    )


def futures_rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
        *,
        dr: float = 0.01
) -> float:
    return (
        price(is_call, S, K, T, r + dr, 0, v) -
        price(is_call, S, K, T, r - dr, 0, v)
    ) / 2

def rho2(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01
) -> float:
    return (
        price(is_call, S, K, T, r, b - db, v) -
        price(is_call, S, K, T, r, b + db, v)
    ) / 2

def strike_delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, K + dS, T, r, b, v) -
        price(is_call, S, K - dS, T, r, b, v)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, K + dS, T, r, b, v) -
        2 * price(is_call, S, K, T, r, b, v) +
        price(is_call, S, K - dS, T, r, b, v)
    ) / dS ** 2
