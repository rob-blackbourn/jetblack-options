# jetblack-options

This repository contains reference implementations of option pricing formulae implemented in Python.

The implementations are covered by tests in an attempt to demonstrate their veracity.

Where possible analytic (or "closed form") solutions are tested against numerical (or "finite difference" algorithms) to establish
the efficacy, accuracy and stability of the algorithms.

## Contributions

Contributions are welcome!

The goals of the project are centred on clarity and accuracy, so optimisations or integrations are not useful.

Valuable contributions include:

* More tests. Note that tests are code too, and have to be maintained. We should aim for the smallest complete set possible. Requests to delete are as relevant as those to add.
* More implementations. Pricing models should come with tests, and a numerical bumping framework (where there are no closed form solutions).

## Usage

This library isn't published as a package, and is not meant to be used directly. Copy the code you want, or use a library that derives from this.

However trying the code is quite simple, no external libraries are required. Simply clone the repository, install, and create a new python file in the "scratch folder". An obvious place to start would be with Black-Scholes.

```python
from jetblack_options.european.black_scholes_merton import (
    price,
    delta
)

is_call = True
S = 110 # Stock price.
K = 100 # Strike price.
r = 0.1 # 10% risk free rate.
q = 0.08 # 8% dividend.
T = 6/12 # Half a year till expiry.
v = 0.125 # 12,5% volatility.

b = r - q # Cost of carry for generalized Black-Scholes.

p = price(is_call, S, K, T, r, b, v)
d = delta(is_call, S, K, T, r, b, v)

# Calculate the delta by bumping the price.
from jetblack_options.numeric_greeks import NumericGreeks
ng = NumericGreeks(price)
d1 = ng.delta(is_call, S, K, T, r, b, v)
```

Typically each module has a function called `price` which calculates the price of the model. Where closed form
solutions are available, they have their canonical name (e.g. `delta`). In this case a general numeric
class was also available.

A tree solution only provides a price. However the numeric class can provide support for the greeks.

```python
from jetblack_options.trees.cox_ross_rubinstein import greeks
from jetblack_options.numeric_greeks import NumericGreeks

is_european = True
is_call = True
S = 110 # Stock price.
K = 100 # Strike price.
r = 0.1 # 10% risk free rate.
q = 0.08 # 8% dividend.
T = 6/12 # Half a year till expiry.
v = 0.125 # 12,5% volatility.

b = r - q
value, delta, gamma, theta = greeks(is_european, is_call, S, K, T, r, b, v, 200)

# Make a lambda to handle `is_european``.
ng = NumericGreeks(lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0])
# The numeric delta should be close to the analytic.
numeric_delta = ng.delta(is_call, S, K, T, r, b, v)
```

## Distributions

Many option pricing formula require probability distribution functions. Sample  implementations are provided,
but where possible that of the standard library `statistics.NormalDist` is used.

Algorithms taking probability functions should provide these as optional arguments to allow the
testing framework to establish the sensitivity to different distribution implementations.
