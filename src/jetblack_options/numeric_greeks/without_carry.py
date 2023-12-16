"""Class for calculating numeric greeks for options using finite difference
methods for the style with no carry or dividend yield.
"""

from typing import Callable, Literal

OptionValue = Callable[
    [
        float,  # Asset price.
        float,  # Strike.
        float,  # Time to expiry in years.
        float,  # Risk free rate.
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
            v: float,
            *,
            dS: float = 0.01,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S + dS, K, T, r, v)
                - self.price(S - dS, K, T, r, v)
            ) / (2 * dS)
        elif method == 'forward':
            return (
                self.price(S + dS, K, T, r, v)
                - self.price(S, K, T, r, v)
            ) / dS
        elif method == 'backward':
            return (
                self.price(S, K, T, r, v)
                - self.price(S - dS, K, T, r, v)
            ) / dS
        else:
            raise ValueError("Invalid method")

    def gamma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S + dS, K, T, r, v)
                - 2 * self.price(S, K, T, r, v)
                + self.price(S - dS, K, T, r, v)
            ) / dS ** 2
        elif method == 'forward':
            return (
                self.price(S + 2 * dS, K, T, r, v)
                - 2 * self.price(S + dS, K, T, r, v)
                + self.price(S, K, T, r, v)
            ) / dS ** 2
        if method == 'central':
            return (
                self.price(S, K, T, r, v)
                - 2 * self.price(S - dS, K, T, r, v)
                + self.price(S - 2 * dS, K, T, r, v)
            ) / dS ** 2
        else:
            raise ValueError("Invalid method")

    def theta(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dT: float = 1 / 365,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T - dT, r, v)
                - self.price(S, K, T + dT, r, v)
            ) / (2 * dT)
        elif method == 'forward':
            return (
                self.price(S, K, T, r, v)
                - self.price(S, K, T + dT, r, v)
            ) / dT
        elif method == 'backward':
            return (
                self.price(S, K, T - dT, r, v)
                - self.price(S, K, T, r, v)
            ) / dT
        else:
            raise ValueError("Invalid method")

    def vega(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dv: float = 0.001,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T, r, v + dv)
                - self.price(S, K, T, r, v - dv)
            ) / (2 * dv)
        elif method == 'forward':
            return (
                self.price(S, K, T, r, v + dv)
                - self.price(S, K, T, r, v)
            ) / dv
        elif method == 'backward':
            return (
                self.price(S, K, T, r, v)
                - self.price(S, K, T, r, v - dv)
            ) / dv
        else:
            raise ValueError('Invalid method')

    def rho(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dr: float = 0.001,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T, r + dr, v)
                - self.price(S, K, T, r - dr, v)
            ) / (2 * dr)
        elif method == 'forward':
            return (
                self.price(S, K, T, r + dr, v)
                - self.price(S, K, T, r, v)
            ) / dr
        elif method == 'backward':
            return (
                self.price(S, K, T, r, v)
                - self.price(S, K, T, r - dr, v)
            ) / dr
        else:
            raise ValueError('Invalid method')

    def elasticity(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
            method: DifferenceMethod = 'central'
    ) -> float:
        return (
            self.delta(S, K, T, r, v, dS=dS, method=method)
            * S
            / self.price(S, K, T, r, v)
        )

    def speed(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(S + 2 * dS, K, T, r, v)
            - 3 * self.price(S + dS, K, T, r, v)
            + 3 * self.price(S, K, T, r, v)
            - self.price(S - dS, K, T, r, v)
        ) / dS ** 3

    def deltap(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(S * (1 + dS), K, T, r, v)
            - self.price(S * (1 - dS), K, T, r, v)
        ) / (2 * S)

    def gammap(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
            method: DifferenceMethod = 'central'
    ) -> float:
        return S / 100 * self.gamma(S, K, T, r, v, dS=dS, method=method)

    def vegap(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(S, K, T, r, v + dv)
            - self.price(S, K, T, r, v - dv)
        ) * v / 0.1 / 2

    def vanna(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
            dv: float = 0.001
    ) -> float:
        # Also known as DdeltaDvol
        return (
            self.price(S + dS, K, T, r, v + dv)
            - self.price(S + dS, K, T, r, v - dv)
            - self.price(S - dS, K, T, r, v + dv)
            + self.price(S - dS, K, T, r, v - dv)
        ) / (4 * dS) / dv

    def charm(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
            dT: float = 1 / 365
    ) -> float:
        # Also known as DdeltaDtime
        return (
            self.price(S + dS, K, T + dT, r, v)
            - self.price(S + dS, K, T - dT, r, v)
            - self.price(S - dS, K, T + dT, r, v)
            + self.price(S - dS, K, T - dT, r, v)
        ) / (4 * dS) / -dT

    def dgamma_dvol(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
            dv: float = 0.01
    ) -> float:
        return (
            self.price(S + dS, K, T, r, v + dv)
            - 2 * self.price(S, K, T, r, v + dv)
            + self.price(S - dS, K, T, r, v + dv)
            - self.price(S + dS, K, T, r, v - dv)
            + 2 * self.price(S, K, T, r, v - dv)
            - self.price(S - dS, K, T, r, v - dv)
        ) / (2 * dv * dS ** 2) / 100

    def vomma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dv: float = 0.001,
    ) -> float:
        # DvegaDvol
        return (
            self.price(S, K, T, r, v + dv)
            - 2 * self.price(S, K, T, r, v)
            + self.price(S, K, T, r, v - dv)
        ) / dv ** 2

    def time_gamma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dT: float = 1 / 365,
    ) -> float:
        return (
            self.price(S, K, T + dT, r, v)
            - 2 * self.price(S, K, T, r, v)
            + self.price(S, K, T - dT, r, v)
        ) / dT ** 2

    def strike_delta(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dX: float = 0.01,
    ) -> float:
        return (
            self.price(S, K + dX, T, r, v)
            - self.price(S, K - dX, T, r, v)
        ) / (2 * dX)

    def strike_gamma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dX: float = 0.01,
    ) -> float:
        return (
            self.price(S, K + dX, T, r, v)
            - 2 * self.price(S, K, T, r, v)
            + self.price(S, K - dX, T, r, v)
        ) / dX ** 2
