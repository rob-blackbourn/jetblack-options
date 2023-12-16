r"""Black (1976) Options on futures/forwards


* The discounted futures price $ F $,
* Strike price $ K $,
* Risk-free rate $ r $,
* Annual dividend yield $ q $,
* Time to maturity $ \tau = T - t $
* Volatility $ \sigma $.

Most of the formula use one or both of the following terms.

$$
d_1 = \frac{\ln(F/K) + (\sigma^2/2)T}{\sigma\sqrt{T}}
$$

$$
d_2 = \frac{\ln(F/K) - (\sigma^2/2)T}{\sigma\sqrt{T}} = d_1 - \sigma\sqrt{T}
$$

$$
\varphi(x) &= \frac{1}{\sqrt{2\pi}} e^{-\frac{1}{2} x^2}
$$

$$
\Phi(x) &= \frac{1}{\sqrt{2\pi}} \int_{-\infty}^x e^{-\frac{1}{2} y^2} \,dy = 1 - \frac{1}{\sqrt{2\pi}} \int_x^\infty e^{-\frac{1}{2} y^2} \,dy
$$
"""

from math import exp, log, sqrt
from statistics import NormalDist

from ..implied_volatility import solve_ivol
from ..numeric_greeks.without_carry import NumericGreeks

norm = NormalDist()
cdf = norm.cdf
pdf = norm.pdf


def price(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    r"""Fair value of a futures/forward using Black 76.

    Call price: $ e^{-r \tau}[F\Phi(d_1) - K\Phi(d_2)] $

    Put price: $ e^{-r \tau} [K\Phi(-d_2) -  F\Phi(-d_1)] $

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The price of the future.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.

    Returns:
        float: The option price.
    """

    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return exp(-r * T) * (F * cdf(d1) - K * cdf(d2))
    else:
        return exp(-r * T) * (K * cdf(-d2) - F * cdf(-d1))


def ivol(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        p: float,
        *,
        max_iterations: int = 20,
        epsilon=1e-8
) -> float:
    r"""Calculate the volatility of a Black 76 option that is implied by the price.

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        p (float): The option price.
        max_iterations (int, Optional): The maximum number of iterations before
            a price is returned. Defaults to 20.
        epsilon (float, Optional): The largest acceptable error. Defaults to 1e-8.

    Returns:
        float: The implied volatility.
    """
    return solve_ivol(
        p,
        lambda v: price(is_call, F, K, T, r, v),
        max_iterations=max_iterations,
        epsilon=epsilon
    )


def make_numeric_greeks(is_call: bool) -> NumericGreeks:
    def evaluate(S: float, K: float, T: float, r: float, b: float) -> float:
        return price(is_call, S, K, T, r, b)
    return NumericGreeks(evaluate)


def delta(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    r"""The sensitivity of the option to a change in the asset price
    using Black 76.

    Call $\Delta$: $e^{-r \tau} \Phi(d_1)$

    Put $\Delta$: $-e^{-r \tau} \Phi(-d_1)$

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The current futures price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The volatility.

    Returns:
        float: The delta.
    """
    d1 = (log(F / K) + T * (v ** 2 / 2)) / (v * sqrt(T))
    if is_call:
        return exp(-r * T) * cdf(d1)
    else:
        return -exp(-r * T) * cdf(-d1)


def gamma(
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    r"""The second derivative to the change in asset price using Black 76.

    $$
    e^{-r \tau} \frac{\varphi(d_1)}{F\sigma\sqrt{\tau}} = K e^{-r \tau} \frac{\varphi(d_2)}{F^2\sigma\sqrt{\tau}}
    $$

    Args:
        F (float): The current futures price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The volatility.

    Returns:
        float: The gamma.
    """
    d1 = (log(F / K) + T * (v ** 2 / 2)) / (v * sqrt(T))

    return exp(-r * T) * pdf(d1) / (F * v * sqrt(T))


def theta(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    r"""The change in the value of the option with respect to time to expiry
    using Black 76.

    $$
    - \frac{F e^{-r \tau} \varphi(d_1) \sigma}{2 \sqrt{\tau}} - rKe^{-r \tau}\Phi(d_2) + rFe^{-r \tau}\Phi(d_1)
    $$

    $$
    - \frac{F e^{-r \tau} \varphi(d_1) \sigma}{2 \sqrt{\tau}} + rKe^{-r \tau}\Phi(-d_2) - rFe^{-r \tau}\Phi(-d_1)
    $$

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The current futures price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The volatility.

    Returns:
        float: The theta.
    """
    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return (
            -F * exp(-r * T) * pdf(d1) * v / (2 * sqrt(T))
            + r * F * exp(-r * T) * cdf(d1)
            - r * K * exp(-r * T) * cdf(d2)
        )
    else:
        return (
            -F * exp(-r * T) * pdf(d1) * v / (2 * sqrt(T))
            - r * F * exp(-r * T) * cdf(-d1)
            + r * K * exp(-r * T) * cdf(-d2)
        )


def vega(
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    r"""The sensitivity of the options price or a change in the asset volatility
    using Black 76.

    $$
    F e^{-r \tau} \varphi(d_1) \sqrt{\tau} = K e^{-r \tau} \varphi(d_2) \sqrt{\tau}
    $$

    Args:
        F (float): The current futures price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The volatility.

    Returns:
        float: The vega.
    """
    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    return F * exp(-r * T) * pdf(d1) * sqrt(T)


def rho(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    r"""The sensitivity of the option price to a change in the risk free rate
    using Black 76.

    Call $ \rho = -\tau e^{-r \tau}[F\Phi(d_1) - K\Phi(d_2)] $

    Put $ \rho = -\tau e^{-r \tau} [K\Phi(-d_2) -  F\Phi(-d_1)] $

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The price of the future.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.

    Returns:
        float: The rho.
    """
    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return -T * exp(-r * T) * (F * cdf(d1) - K * cdf(d2))
    else:
        return -T * exp(-r * T) * (K * cdf(-d2) - F * cdf(-d1))


def vanna(
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    r"""The sensitivity of the option value to the underlying
    asset price and the volatility.

    $$
    \frac{\partial^2 V}{\partial F \partial \sigma} = -e^{-r \tau} \varphi(d_1) \frac{d_2}{\sigma} \, = \frac{\mathcal{V}}{F}\left[1 - \frac{d_1}{\sigma\sqrt{\tau}} \right]
    $$

    Args:
        F (float): The price of the future.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.

    Returns:
        float: The vanna.
    """
    d1 = (log(F / K) + T * (v ** 2 / 2)) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return -exp(-r * T) * pdf(d1) * d2 / v


def vomma(
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
) -> float:
    r"""The second order sensitivity to volatility.

    $$
    F e^{-r \tau} \varphi(d_1) \sqrt{\tau} \frac{d_1 d_2}{\sigma} = \mathcal{V}  \frac{d_1 d_2}{\sigma}
    $$

    Args:
        F (float): The price of the future.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.

    Returns:
        float: The vomma
    """
    d1 = (log(F / K) + T * (v ** 2 / 2)) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return F * exp(-r * T) * pdf(d1) * sqrt(T) * d1 * d2 / v
