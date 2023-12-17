An obvious place to start would be with Black, Scholes & Merton. This was the
classic European option pricing formula that provides a method for pricing
European stock options with a dividend yield.

```python
from jetblack_options.european.black_scholes_merton import (
    price,
    delta
)

is_call = True
S = 110 # Asset price.
K = 100 # Strike price.
r = 0.1 # 10% risk free rate.
q = 0.08 # 8% dividend.
T = 6/12 # Half a year till expiry.
v = 0.125 # 12.5% volatility.

p = price(is_call, S, K, T, r, q, v)
d = delta(is_call, S, K, T, r, q, v)
```

## Module contents

Almost all modules will contain a `price` function.

Given a function to calculate the price an `ivol` function will exists to
calculate the implied volatility. This is done with a solver.

If a `NumericalGreeks` class has been written for the parameters used in the
model a `make_numeric_greeks` function will be available
(see [Numeric Greeks](./numeric-greeks.md)).

Finally if closed form solutions are available, functions may have been
written for these. For example a large number of greeks are available for
much of the Black Scholes universe.

## What next ?

[Numeric greeks](./numeric-greeks.md)