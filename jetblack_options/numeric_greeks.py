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

    def dgamma_dvol(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S + dS, X, T, r, b, v + 0.01) -
            2 * self.price(is_call, S, X, T, r, b, v + 0.01) +
            self.price(is_call, S - dS, X, T, r, b, v + 0.01) -
            self.price(is_call, S + dS, X, T, r, b, v - 0.01) +
            2 * self.price(is_call, S, X, T, r, b, v - 0.01) -
            self.price(is_call, S - dS, X, T, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100

    def gammap(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return S / 100 * (
            self.price(is_call, S + dS, X, T, r, b, v) -
            2 * self.price(is_call, S, X, T, r, b, v) +
            self.price(is_call, S - dS, X, T, r, b, v)
        ) / dS ** 2

    def time_gamma(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dT: float = 1 / 365,
    ) -> float:
        return (
            self.price(is_call, S, X, T + dT, r, b, v) -
            2 * self.price(is_call, S, X, T, r, b, v) +
            self.price(is_call, S, X, T - dT, r, b, v)
        ) / dT ** 2

    def ddelta_dvol(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return 1 / (4 * dS * 0.01) * (
            self.price(is_call, S + dS, X, T, r, b, v + 0.01) -
            self.price(is_call, S + dS, X, T, r, b, v - 0.01) -
            self.price(is_call, S - dS, X, T, r, b, v + 0.01) +
            self.price(is_call, S - dS, X, T, r, b, v - 0.01)
        ) / 100

    def vega(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X, T, r, b, v + dv) -
            self.price(is_call, S, X, T, r, b, v - dv)
        ) / 2

    def vomma(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X, T, r, b, v + dv) -
            2 * self.price(is_call, S, X, T, r, b, v) +
            self.price(is_call, S, X, T, r, b, v - dv)
        ) / dv ** 2 / 10000

    def vegap(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X, T, r, b, v + dv) -
            self.price(is_call, S, X, T, r, b, v - dv)
        ) * v / 0.1 / 2

    def dvega_dvol(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X, T, r, b, v + dv) -
            2 * self.price(is_call, S, X, T, r, b, v) +
            self.price(is_call, S, X, T, r, b, v - dv)
        )

    def theta(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dT: float = 1 /365,
    ) -> float:
        if T <= dT:
            return (
                self.price(is_call, S, X, 0.00001, r, b, v)
                - self.price(is_call, S, X, T, r, b, v)
            )
        else:
            return (
                self.price(is_call, S, X, T - dT, r, b, v)
                - self.price(is_call, S, X, T, r, b, v)
            )

    def rho(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dr: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X, T, r + dr, b + dr, v)
            - self.price(is_call, S, X, T, r - dr, b - dr, v)
        ) / 2

    def futures_rho(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dr: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X, T, r + dr, b, v)
            - self.price(is_call, S, X, T, r - dr, b, v)
        ) / 2

    def rho2(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            db: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X, T, r, b - db, v)
            - self.price(is_call, S, X, T, r, b + db, v)
        ) / 2

    def carry(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            db: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X, T, r, b + db, v)
            - self.price(is_call, S, X, T, r, b - db, v)
        ) / 2

    def speed(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S + 2 * dS, X, T, r, b, v)
            - 3 * self.price(is_call, S + dS, X, T, r, b, v)
            + 3 * self.price(is_call, S, X, T, r, b, v)
            - self.price(is_call, S - dS, X, T, r, b, v)
        ) / dS ** 3

    def strike_delta(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X + dS, T, r, b, v)
            - self.price(is_call, S, X - dS, T, r, b, v)
        ) / (2 * dS)

    def strike_gamma(
            self,
            is_call: bool,
            S: float,
            X: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return (
            self.price(is_call, S, X + dS, T, r, b, v)
            - 2 * self.price(is_call, S, X, T, r, b, v)
            + self.price(is_call, S, X - dS, T, r, b, v)
        ) / dS ** 2
