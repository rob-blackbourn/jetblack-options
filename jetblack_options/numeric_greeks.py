"""Class for calculating numeric greeks for an option"""

from typing import Callable

OptionValue = Callable[
    [
        bool, # is_call
        float, # Asset price.
        float, # Strike.
        float, # Time to expiry in years.
        float, # Risk free rate.
        float, # Cost of carry.
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
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
                self.price(is_call, S + dS, K, T, r, b, v) -
                self.price(is_call, S - dS, K, T, r, b, v)
        ) / (2 * dS)

    def gamma(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S + dS, K, T, r, b, v) -
            2 * self.price(is_call, S, K, T, r, b, v) +
            self.price(is_call, S - dS, K, T, r, b, v)
        ) / dS ** 2

    def theta(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dT: float = 1 /365,
    ) -> float:
        if T <= dT:
            return (
                self.price(is_call, S, K, 0.00001, r, b, v)
                - self.price(is_call, S, K, T, r, b, v)
            )
        else:
            return (
                self.price(is_call, S, K, T - dT, r, b, v)
                - self.price(is_call, S, K, T, r, b, v)
            )

    def vega(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r, b, v + dv) -
            self.price(is_call, S, K, T, r, b, v - dv)
        ) / 2

    def rho(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dr: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r + dr, b + dr, v)
            - self.price(is_call, S, K, T, r - dr, b - dr, v)
        ) / 2

    def carry(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            db: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r, b + db, v)
            - self.price(is_call, S, K, T, r, b - db, v)
        ) / 2

    def elasticity(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S + dS, K, T, r, b, v) -
            self.price(is_call, S - dS, K, T, r, b, v)
        ) / (2 * dS) * S / self.price(is_call, S, K, T, r, b, v)

    def speed(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S + 2 * dS, K, T, r, b, v)
            - 3 * self.price(is_call, S + dS, K, T, r, b, v)
            + 3 * self.price(is_call, S, K, T, r, b, v)
            - self.price(is_call, S - dS, K, T, r, b, v)
        ) / dS ** 3


    def deltap(
            self,
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
            self.price(is_call, S * (1 + dS), K, T, r, b, v) -
            self.price(is_call, S * (1 - dS), K, T, r, b, v)
        ) * 2 / S

    def gammap(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return S / 100 * (
            self.price(is_call, S + dS, K, T, r, b, v) -
            2 * self.price(is_call, S, K, T, r, b, v) +
            self.price(is_call, S - dS, K, T, r, b, v)
        ) / dS ** 2

    def vegap(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r, b, v + dv) -
            self.price(is_call, S, K, T, r, b, v - dv)
        ) * v / 0.1 / 2

    def ddelta_dvol(
            self,
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
        # Also known as vanna
        return 1 / (4 * dS * 0.01) * (
            self.price(is_call, S + dS, K, T, r, b, v + dv) -
            self.price(is_call, S + dS, K, T, r, b, v - dv) -
            self.price(is_call, S - dS, K, T, r, b, v + dv) +
            self.price(is_call, S - dS, K, T, r, b, v - dv)
        ) / 100

    def dgamma_dvol(
            self,
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
            self.price(is_call, S + dS, K, T, r, b, v + dv) -
            2 * self.price(is_call, S, K, T, r, b, v + dv) +
            self.price(is_call, S - dS, K, T, r, b, v + dv) -
            self.price(is_call, S + dS, K, T, r, b, v - dv) +
            2 * self.price(is_call, S, K, T, r, b, v - dv) -
            self.price(is_call, S - dS, K, T, r, b, v - dv)
        ) / (2 * dv * dS ** 2) / 100

    def dvega_dvol(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r, b, v + dv) -
            2 * self.price(is_call, S, K, T, r, b, v) +
            self.price(is_call, S, K, T, r, b, v - dv)
        )

    def vomma(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
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
            b,
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
            b: float,
            v: float,
            *,
            dT: float = 1 / 365,
    ) -> float:
        return (
            self.price(is_call, S, K, T + dT, r, b, v) -
            2 * self.price(is_call, S, K, T, r, b, v) +
            self.price(is_call, S, K, T - dT, r, b, v)
        ) / dT ** 2

    def futures_rho(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dr: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r + dr, b, v)
            - self.price(is_call, S, K, T, r - dr, b, v)
        ) / 2

    def rho2(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            db: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K, T, r, b - db, v)
            - self.price(is_call, S, K, T, r, b + db, v)
        ) / 2

    def strike_delta(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dX: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K + dX, T, r, b, v)
            - self.price(is_call, S, K - dX, T, r, b, v)
        ) / (2 * dX)

    def strike_gamma(
            self,
            is_call: bool,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dX: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, K + dX, T, r, b, v)
            - 2 * self.price(is_call, S, K, T, r, b, v)
            + self.price(is_call, S, K - dX, T, r, b, v)
        ) / dX ** 2
