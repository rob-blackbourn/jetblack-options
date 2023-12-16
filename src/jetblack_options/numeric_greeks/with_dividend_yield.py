"""Class for calculating numeric greeks for options using finite difference
methods for the dividend yield style.
"""

from typing import Callable, Literal

OptionValue = Callable[
    [
        float,  # Asset price.
        float,  # Strike.
        float,  # Time to expiry in years.
        float,  # Risk free rate.
        float,  # Dividend yield.
        float  # Asset volatility
    ],
    float  # The option price
]
DifferenceMethod = Literal['central', 'forward', 'backward']


class NumericGreeks:

    def __init__(
            self,
            price: OptionValue
    ) -> None:
        self.price = price

    def delta(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S + dS, K, T, r, q, v)
                - self.price(S - dS, K, T, r, q, v)
            ) / (2 * dS)
        elif method == 'forward':
            return (
                self.price(S + dS, K, T, r, q, v)
                - self.price(S, K, T, r, q, v)
            ) / dS
        elif method == 'backward':
            return (
                self.price(S, K, T, r, q, v)
                - self.price(S - dS, K, T, r, q, v)
            ) / dS
        else:
            raise ValueError("Invalid method")

    def gamma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S + dS, K, T, r, q, v)
                - 2 * self.price(S, K, T, r, q, v)
                + self.price(S - dS, K, T, r, q, v)
            ) / dS ** 2
        elif method == 'central':
            return (
                self.price(S + 2 * dS, K, T, r, q, v)
                - 2 * self.price(S + dS, K, T, r, q, v)
                + self.price(S, K, T, r, q, v)
            ) / dS ** 2
        if method == 'backward':
            return (
                self.price(S, K, T, r, q, v)
                - 2 * self.price(S - dS, K, T, r, q, v)
                + self.price(S - 2 * dS, K, T, r, q, v)
            ) / dS ** 2
        else:
            raise ValueError("Invalid method")

    def theta(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dT: float = 1 / 365,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T - dT, r, q, v)
                - self.price(S, K, T + dT, r, q, v)
            ) / (2 * dT)
        elif method == 'forward':
            return (
                self.price(S, K, T, r, q, v)
                - self.price(S, K, T + dT, r, q, v)
            ) / dT
        elif method == 'backward':
            return (
                self.price(S, K, T - dT, r, q, v)
                - self.price(S, K, T, r, q, v)
            ) / dT
        else:
            raise ValueError("Invalid method")

    def vega(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dv: float = 0.001,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T, r, q, v + dv)
                - self.price(S, K, T, r, q, v - dv)
            ) / (2 * dv)
        elif method == 'forward':
            return (
                self.price(S, K, T, r, q, v + dv)
                - self.price(S, K, T, r, q, v)
            ) / dv
        elif method == 'backward':
            return (
                self.price(S, K, T, r, q, v)
                - self.price(S, K, T, r, q, v - dv)
            ) / dv
        else:
            raise ValueError('Invalid method')

    def rho(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dr: float = 0.001,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T, r + dr, q, v)
                - self.price(S, K, T, r - dr, q, v)
            ) / (2 * dr)
        elif method == 'forward':
            return (
                self.price(S, K, T, r + dr, q, v)
                - self.price(S, K, T, r, q, v)
            ) / dr
        elif method == 'backward':
            return (
                self.price(S, K, T, r, q, v)
                - self.price(S, K, T, r - dr, q, v)
            ) / dr
        else:
            raise ValueError('Invalid method')

    def carry(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dq: float = 0.001,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T, r, q + dq, v)
                - self.price(S, K, T, r, q - dq, v)
            ) / (2 * dq)
        elif method == 'forward':
            return (
                self.price(S, K, T, r, q + dq, v)
                - self.price(S, K, T, r, q, v)
            ) / dq
        elif method == 'backward':
            return (
                self.price(S, K, T, r, q, v)
                - self.price(S, K, T, r, q - dq, v)
            ) / dq
        else:
            raise ValueError('Invalid method')

    def elasticity(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(S + dS, K, T, r, q, v)
            - self.price(S - dS, K, T, r, q, v)
        ) / (2 * dS) * S / self.price(S, K, T, r, q, v)

    def speed(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(S + 2 * dS, K, T, r, q, v)
            - 3 * self.price(S + dS, K, T, r, q, v)
            + 3 * self.price(S, K, T, r, q, v)
            - self.price(S - dS, K, T, r, q, v)
        ) / dS ** 3

    def deltap(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01
    ) -> float:
        return (
            self.price(S * (1 + dS), K, T, r, q, v)
            - self.price(S * (1 - dS), K, T, r, q, v)
        ) * 2 / S

    def gammap(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return S / 100 * (
            self.price(S + dS, K, T, r, q, v)
            - 2 * self.price(S, K, T, r, q, v)
            + self.price(S - dS, K, T, r, q, v)
        ) / dS ** 2

    def vegap(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(S, K, T, r, q, v + dv)
            - self.price(S, K, T, r, q, v - dv)
        ) * v / 0.1 / 2

    def vanna(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01,
            dv: float = 0.001
    ) -> float:
        # Also known as DdeltaDvol
        return (
            self.price(S + dS, K, T, r, q, v + dv)
            - self.price(S + dS, K, T, r, q, v - dv)
            - self.price(S - dS, K, T, r, q, v + dv)
            + self.price(S - dS, K, T, r, q, v - dv)
        ) / (4 * dS) / dv

    def charm(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01,
            dT: float = 1 / 365
    ) -> float:
        # Also known as DdeltaDtime
        return (
            self.price(S + dS, K, T + dT, r, q, v)
            - self.price(S + dS, K, T - dT, r, q, v)
            - self.price(S - dS, K, T + dT, r, q, v)
            + self.price(S - dS, K, T - dT, r, q, v)
        ) / (4 * dS) / -dT

    def dgamma_dvol(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dS: float = 0.01,
            dv: float = 0.001
    ) -> float:
        return (
            self.price(S + dS, K, T, r, q, v + dv)
            - 2 * self.price(S, K, T, r, q, v + dv)
            + self.price(S - dS, K, T, r, q, v + dv)
            - self.price(S + dS, K, T, r, q, v - dv)
            + 2 * self.price(S, K, T, r, q, v - dv)
            - self.price(S - dS, K, T, r, q, v - dv)
        ) / (2 * dv * dS ** 2)

    def vomma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dv: float = 0.001,
    ) -> float:
        # DvegaDvol
        return (
            self.price(S, K, T, r, q, v + dv)
            - 2 * self.price(S, K, T, r, q, v)
            + self.price(S, K, T, r, q, v - dv)
        ) / dv ** 2

    def time_gamma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dT: float = 1 / 365,
    ) -> float:
        return (
            self.price(S, K, T + dT, r, q, v)
            - 2 * self.price(S, K, T, r, q, v)
            + self.price(S, K, T - dT, r, q, v)
        ) / dT ** 2

    def futures_rho(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dr: float = 0.01,
    ) -> float:
        return (
            self.price(S, K, T, r + dr, q, v)
            - self.price(S, K, T, r - dr, q, v)
        ) / 2

    def rho2(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dq: float = 0.01,
    ) -> float:
        return (
            self.price(S, K, T, r, q - dq, v)
            - self.price(S, K, T, r, q + dq, v)
        ) / 2

    def strike_delta(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dK: float = 0.01,
    ) -> float:
        return (
            self.price(S, K + dK, T, r, q, v)
            - self.price(S, K - dK, T, r, q, v)
        ) / (2 * dK)

    def strike_gamma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            q: float,
            v: float,
            *,
            dK: float = 0.01,
    ) -> float:
        return (
            self.price(S, K + dK, T, r, q, v)
            - 2 * self.price(S, K, T, r, q, v)
            + self.price(S, K - dK, T, r, q, v)
        ) / dK ** 2
