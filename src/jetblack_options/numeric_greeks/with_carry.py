"""Class for calculating numeric greeks for options using finite difference
methods for the generalised style using cost of carry.
"""

from typing import Callable, Literal

OptionValue = Callable[
    [
        float,  # Asset price.
        float,  # Strike.
        float,  # Time to expiry in years.
        float,  # Risk free rate.
        float,  # Cost of carry.
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
            b: float,
            v: float,
            *,
            dS: float = 0.01,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S + dS, K, T, r, b, v)
                - self.price(S - dS, K, T, r, b, v)
            ) / (2 * dS)
        elif method == 'forward':
            return (
                self.price(S + dS, K, T, r, b, v)
                - self.price(S, K, T, r, b, v)
            ) / dS
        elif method == 'backward':
            return (
                self.price(S, K, T, r, b, v)
                - self.price(S - dS, K, T, r, b, v)
            ) / dS
        else:
            raise ValueError("Invalid method")

    def gamma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S + dS, K, T, r, b, v)
                - 2 * self.price(S, K, T, r, b, v)
                + self.price(S - dS, K, T, r, b, v)
            ) / dS ** 2
        elif method == 'forward':
            return (
                self.price(S + 2 * dS, K, T, r, b, v)
                - 2 * self.price(S + dS, K, T, r, b, v)
                + self.price(S, K, T, r, b, v)
            ) / dS ** 2
        elif method == 'backward':
            return (
                self.price(S, K, T, r, b, v)
                - 2 * self.price(S - dS, K, T, r, b, v)
                + self.price(S - 2 * dS, K, T, r, b, v)
            ) / dS ** 2
        else:
            raise ValueError("Invalid method")

    def theta(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dT: float = 1 / 365,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T - dT, r, b, v)
                - self.price(S, K, T + dT, r, b, v)
            ) / (2 * dT)
        elif method == 'forward':
            return (
                self.price(S, K, T, r, b, v)
                - self.price(S, K, T + dT, r, b, v)
            ) / dT
        if method == 'backward':
            return (
                self.price(S, K, T - dT, r, b, v)
                - self.price(S, K, T, r, b, v)
            ) / dT
        else:
            raise ValueError("Invalid method")

    def vega(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.001,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T, r, b, v + dv)
                - self.price(S, K, T, r, b, v - dv)
            ) / (2 * dv)
        elif method == 'forward':
            return (
                self.price(S, K, T, r, b, v + dv)
                - self.price(S, K, T, r, b, v)
            ) / dv
        elif method == 'backward':
            return (
                self.price(S, K, T, r, b, v)
                - self.price(S, K, T, r, b, v - dv)
            ) / dv
        else:
            raise ValueError('Invalid method')

    def rho(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dr: float = 0.001,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T, r + dr, b + dr, v)
                - self.price(S, K, T, r - dr, b - dr, v)
            ) / (2 * dr)
        elif method == 'central':
            return (
                self.price(S, K, T, r + dr, b + dr, v)
                - self.price(S, K, T, r - dr, b, v)
            ) / dr
        if method == 'backward':
            return (
                self.price(S, K, T, r + dr, b, v)
                - self.price(S, K, T, r - dr, b - dr, v)
            ) / dr
        else:
            raise ValueError('Invalid method')

    def carry(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            db: float = 0.001,
            method: DifferenceMethod = 'central'
    ) -> float:
        if method == 'central':
            return (
                self.price(S, K, T, r, b + db, v)
                - self.price(S, K, T, r, b - db, v)
            ) / (2 * db)
        if method == 'forward':
            return (
                self.price(S, K, T, r, b + db, v)
                - self.price(S, K, T, r, b, v)
            ) / db
        if method == 'central':
            return (
                self.price(S, K, T, r, b, v)
                - self.price(S, K, T, r, b - db, v)
            ) / db
        else:
            raise ValueError('Invalid method')

    def elasticity(
            self,
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
            self.delta(S, K, T, r, b, v, dS=dS) * S
            / self.price(S, K, T, r, b, v)
        )

    def speed(
            self,
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
            self.price(S + 2 * dS, K, T, r, b, v)
            - 3 * self.price(S + dS, K, T, r, b, v)
            + 3 * self.price(S, K, T, r, b, v)
            - self.price(S - dS, K, T, r, b, v)
        ) / dS ** 3

    def deltap(
            self,
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
            self.price(S * (1 + dS), K, T, r, b, v)
            - self.price(S * (1 - dS), K, T, r, b, v)
        ) * 2 / S

    def gammap(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
    ) -> float:
        return self.gamma(S, K, T, r, b, v, dS=dS) * S / 100

    def vegap(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.001,
    ) -> float:
        return self.vega(S, K, T, r, b, v, dv=dv) * v * 10

    def vanna(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
            dv: float = 0.001
    ) -> float:
        """The second order derivative of the option price to a change in the asset
        price and a change in the volatility.

        Args:
            S (float): The asset price.
            K (float): The strike price.
            T (float): The time to expiry in years.
            r (float): The risk free rate.
            b (float): The cost of carry.
            v (float): The asset volatility.
            dS (float, optional): The change in spot price. Defaults to 0.01.
            dv (float, optional): The change in volatility. Defaults to 0.01.

        Returns:
            float: _description_
        """
        # Also known as DdeltaDvol
        return (
            self.price(S + dS, K, T, r, b, v + dv)
            - self.price(S + dS, K, T, r, b, v - dv)
            - self.price(S - dS, K, T, r, b, v + dv)
            + self.price(S - dS, K, T, r, b, v - dv)
        ) / (4 * dS) / dv

    def charm(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
            dT: float = 1 / 365
    ) -> float:
        """Measures the instantaneous rate of change of delta over the passage of
        time.

        Also known as DdeltaDtime.

        Args:
            S (float): The asset price.
            K (float): The strike price.
            T (float): The time to expiry in years.
            r (float): The risk free rate.
            b (float): The cost of carry.
            v (float): The asset volatility.
            dS (float, optional): Change in asset price. Defaults to 0.01.
            dT (float, optional): Change in time. Defaults to 1/365.

        Returns:
            float: The charm.
        """
        # Also known as DdeltaDtime
        return (
            self.price(S + dS, K, T + dT, r, b, v)
            - self.price(S + dS, K, T - dT, r, b, v)
            - self.price(S - dS, K, T + dT, r, b, v)
            + self.price(S - dS, K, T - dT, r, b, v)
        ) / (4 * dS) / -dT

    def dgamma_dvol(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dS: float = 0.01,
            dv: float = 0.001
    ) -> float:
        return (
            self.price(S + dS, K, T, r, b, v + dv)
            - 2 * self.price(S, K, T, r, b, v + dv)
            + self.price(S - dS, K, T, r, b, v + dv)
            - self.price(S + dS, K, T, r, b, v - dv)
            + 2 * self.price(S, K, T, r, b, v - dv)
            - self.price(S - dS, K, T, r, b, v - dv)
        ) / (2 * dv * dS ** 2)

    def vomma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dv: float = 0.001,
    ) -> float:
        # DvegaDvol
        return (
            self.price(S, K, T, r, b, v + dv)
            - 2 * self.price(S, K, T, r, b, v)
            + self.price(S, K, T, r, b, v - dv)
        ) / dv ** 2

    def time_gamma(
            self,
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
            self.price(S, K, T + dT, r, b, v)
            - 2 * self.price(S, K, T, r, b, v)
            + self.price(S, K, T - dT, r, b, v)
        ) / dT ** 2

    def futures_rho(
            self,
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
            self.price(S, K, T, r + dr, b, v)
            - self.price(S, K, T, r - dr, b, v)
        ) / 2

    def rho2(
            self,
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
            self.price(S, K, T, r, b - db, v)
            - self.price(S, K, T, r, b + db, v)
        ) / 2

    def strike_delta(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dK: float = 0.01,
    ) -> float:
        return (
            self.price(S, K + dK, T, r, b, v)
            - self.price(S, K - dK, T, r, b, v)
        ) / (2 * dK)

    def strike_gamma(
            self,
            S: float,
            K: float,
            T: float,
            r: float,
            b: float,
            v: float,
            *,
            dK: float = 0.01,
    ) -> float:
        return (
            self.price(S, K + dK, T, r, b, v)
            - 2 * self.price(S, K, T, r, b, v)
            + self.price(S, K - dK, T, r, b, v)
        ) / dK ** 2
