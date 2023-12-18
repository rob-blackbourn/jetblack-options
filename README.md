# jetblack-options

This repository contains reference implementations of option pricing formulae
implemented in Python.

It has no dependencies.

There is an online demonstration of some of the valuations [here](https://rob-blackbourn.github.io/jetblack-options-demo/)
(source code [here](https://github.com/rob-blackbourn/jetblack-options-demo)).

## Status

This is currently considered alpha.

## Usage

The library can be installed as a package.

```bash
pip install jetblack-options
```

An obvious place to start would be with Black-Scholes.

```python
from jetblack_options.european.black_scholes_merton import (
    price,
    delta,
    make_numeric_greeks
)

is_call = True
S = 110 # Asset price.
K = 100 # Strike price.
r = 0.1 # 10% risk free rate.
q = 0.08 # 8% dividend.
T = 6/12 # Half a year till expiry.
v = 0.125 # 12.5% volatility.

b = r - q # Cost of carry for generalized Black-Scholes.

p = price(is_call, S, K, T, r, b, v)
d = delta(is_call, S, K, T, r, b, v)

# Calculate the delta by bumping the price.
ng = make_numeric_greeks(is_call=True)
d1 = ng.delta(S, K, T, r, b, v)
```

For more information [read the docs](https://rob-blackbourn.github.io/jetblack-options/).
