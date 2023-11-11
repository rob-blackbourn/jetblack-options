"""Class for calculating numeric greeks for an option"""

from typing import Callable

OptionValue = Callable[
    [
        bool, # is_call
        float, # Asset price.
        float, # Strike.
        float, # Time to expiry in years.
        float, # Risk free rate.
        float # Asset volatility
    ],
    float # The option price
]

class NumericGreeks:

    def __init__(
            self,
            price: OptionValue
    ) -> None:
        self.price = price

    def delta(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
                self.price(is_call, S + dS, K, T, r, v) -
                self.price(is_call, S - dS, K, T, r, v)
        ) / (2 * dS)

    def gamma(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        p1 = self.price(is_call, S + dS, K, T, r, v)
        p2 = self.price(is_call, S, K, T, r, v)
        p3 = self.price(is_call, S - dS, K, T, r, v)
        return (p1 - 2 * p2 + p3) / dS ** 2

    def theta(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dT: float = 1 /365,
    ) -> float:
        if T <= dT:
            return (
                self.price(is_call, S, K, 0.00001, r, v)
                - self.price(is_call, S, K, T, r, v)
            )
        else:
            return (
                self.price(is_call, S, K, T - dT, r, v)
                - self.price(is_call, S, K, T + dT, r, v)
            ) / 2

    def vega(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r, v + dv)
            - self.price(is_call, S, K, T, r, v - dv)
        ) / 2

    def rho(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dr: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r + dr, v)
            - self.price(is_call, S, K, T, r - dr, v)
        ) / 2

    def elasticity(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S + dS, K, T, r, v) -
            self.price(is_call, S - dS, K, T, r, v)
        ) / (2 * dS) * S / self.price(is_call, S, K, T, r, v)

    def speed(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S + 2 * dS, K, T, r, v)
            - 3 * self.price(is_call, S + dS, K, T, r, v)
            + 3 * self.price(is_call, S, K, T, r, v)
            - self.price(is_call, S - dS, K, T, r, v)
        ) / dS ** 3


    def deltap(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01
    ) -> float:
        return (
            self.price(is_call, S * (1 + dS), K, T, r, v) -
            self.price(is_call, S * (1 - dS), K, T, r, v)
        ) * 2 / S

    def gammap(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return S / 100 * (
            self.price(is_call, S + dS, K, T, r, v) -
            2 * self.price(is_call, S, K, T, r, v) +
            self.price(is_call, S - dS, K, T, r, v)
        ) / dS ** 2

    def vegap(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r, v + dv) -
            self.price(is_call, S, K, T, r, v - dv)
        ) * v / 0.1 / 2

    def vanna(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dS: float = 0.01,
            dv: float = 0.01
    ) -> float:
        # Also known as DdeltaDvol
        return 1 / (4 * dS * 0.01) * (
            self.price(is_call, S + dS, K, T, r, v + dv) -
            self.price(is_call, S + dS, K, T, r, v - dv) -
            self.price(is_call, S - dS, K, T, r, v + dv) +
            self.price(is_call, S - dS, K, T, r, v - dv)
        ) / 100

    def charm(
            self,
            is_call: bool,
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
            self.price(is_call, S + dS, K, T + dT, r, v) -
            self.price(is_call, S + dS, K, T - dT, r, v) -
            self.price(is_call, S - dS, K, T + dT, r, v) +
            self.price(is_call, S - dS, K, T - dT, r, v)
        ) / (4 * dS) / -dT

    def dgamma_dvol(
            self,
            is_call: bool,
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
            self.price(is_call, S + dS, K, T, r, v + dv) -
            2 * self.price(is_call, S, K, T, r, v + dv) +
            self.price(is_call, S - dS, K, T, r, v + dv) -
            self.price(is_call, S + dS, K, T, r, v - dv) +
            2 * self.price(is_call, S, K, T, r, v - dv) -
            self.price(is_call, S - dS, K, T, r, v - dv)
        ) / (2 * dv * dS ** 2) / 100

    def dvega_dvol(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r, v + dv) -
            2 * self.price(is_call, S, K, T, r, v) +
            self.price(is_call, S, K, T, r, v - dv)
        )

    def vomma(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return self.dvega_dvol(
            is_call,
            S,
            K,
            T,
            r,
            v,
            dv=dv
        ) / dv ** 2 / 10000

    def time_gamma(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dT: float = 1 / 365,
    ) -> float:
        return (
            self.price(is_call, S, K, T + dT, r, v) -
            2 * self.price(is_call, S, K, T, r, v) +
            self.price(is_call, S, K, T - dT, r, v)
        ) / dT ** 2

    def strike_delta(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dX: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K + dX, T, r, v)
            - self.price(is_call, S, K - dX, T, r, v)
        ) / (2 * dX)

    def strike_gamma(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            v: float,
            *,
            dX: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K + dX, T, r, v)
            - 2 * self.price(is_call, S, K, T, r, v)
            + self.price(is_call, S, K - dX, T, r, v)
        ) / dX ** 2
