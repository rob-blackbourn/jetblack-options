The "greeks" are measures of the sensitivity of an attribute (often the price)
to one or more of the other attributes. So the *delta* if the sensitivity of
the price of the option to a change in the asset price.

Some option models provide analytical solutions for this. In particular the
Black-Scholes style models have closed form solutions. Many other models, in
particular tree based formulations, do not.

## Finite Difference Methods

The greeks can be calculated numerically by using finite difference methods.
This means calculating the price of the option multiple times, while perturbing
the inputs.

Using the central difference the formula for calculating the delta is given
below.

$$
\frac{\partial C}{\partial S} = \frac{BS_{Call}(S+\Delta S, K, T,r,\sigma) - BS_{Call}(S-\Delta, K, T,r,\sigma)}{2 \Delta S}
$$

The following is a simplistic python implementation.

```python
def delta(is_call, S, K, T, r, v, dS):
    return (
        price(is_call, S + dS, K, T, r, v)
        - price(is_call, S - dS, K, T, r, v)
    ) / (2 * dS)
```

## The `NumericGreeks` classes

Some classes are provided for calculating the greeks. They differ only in the
way they handle carry/dividend yield, or the lack of it.

* `jetblack_options.numeric_greeks.without_carry` - for pricing formula with no
    carry or dividend yield, for example Black 76, or the original Black-Scholes
    formula for non-dividend paying stock.
* `jetblack_options.numeric_greeks.with_dividend_yield` - for pricing formulae
    with a continuous dividend yield.
* `jetblack_options.numeric_greeks.with_carry` - for pricing formulae with cost
    of carry, in the style of the generalised Black Scholes model.

Each of the option pricing models has a convenience method `make_numeric_greeks`
which will choose the appropriate `NumericGreeks` class.

## Optional Arguments

Some of the methods have an optional `method` parameter. This controls which finite difference is used. This can be
one of: `'central'`, `'forward'` or `'backward'`.

All the methods take as an optional parameter the value of
the *bump* being applied. For example the `delta` method
takes a `dS` argument which has the default value of `0.01`.

## Examples

Here we calculate the delta for the Black, Scholes & Merton model with continuous
dividend yield using the finite difference.

```python
# Calculate the delta by bumping the price.
from jetblack_options.european.black_scholes_merton import make_numeric_greeks
ng = make_numeric_greeks(is_call=True)
d1 = ng.delta(is_call, S, K, T, r, q, v)
```


## What next ?

[Using pandas](./pandas.md)