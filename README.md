# jetblack-options

This repository contains reference implementations of option pricing formulae
implemented in Python.

It has no dependencies.

There is a web UI demonstrating some of the valuations [here](https://rob-blackbourn.github.io/jetblack-options/).
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
    delta
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
from jetblack_options.numeric_greeks.with_carry import NumericGreeks
ng = NumericGreeks(price)
d1 = ng.delta(is_call, S, K, T, r, b, v)
```

Typically each module has a function called `price` which calculates the price
of the model. Where closed form solutions are available, they have their
canonical name (e.g. `delta`). In this case a general numeric class was also
available.

A tree solution only provides a price. However the numeric class can provide
support for the greeks.

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

# Make a lambda to handle is_european.
ng = NumericGreeks(
    lambda is_call, S, K, T, r, b, v: greeks(is_european, is_call, S, K, T, r, b, v, 100)[0]
)
# The numeric delta should be close to the analytic.
numeric_delta = ng.delta(is_call, S, K, T, r, b, v)
```

## Pandas

The code has been written without dependencies to keep the implementation clean.

It is fairly simple to support vectors. The following is the generalised
Black Scholes price function using scipy, numpy and pandas.

```python
import numpy as np
import pandas as pd
from scipy.stats import norm

def price(is_call, S, K, T, r, b, v):
    d1 = (np.log(S / K) + T * (b + v ** 2 / 2)) / (v * np.sqrt(T))
    d2 = d1 - v * np.sqrt(T)

    if is_call:
        return S * np.exp((b - r) * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    else:
        return K * np.exp(-r * T) * norm.cdf(-d2) - S * np.exp((b - r) * T) * norm.cdf(-d1)

df = pd.DataFrame([
    {'S': 110, 'K': 100, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'S': 100, 'K': 100, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'S': 100, 'K': 110, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
])

price(True, df['S'], df['K'], df['T'], df['r'], df['r'] - df['q'], df['v'])
```

Here's another implementation using a single data frame.

```python
import numpy as np
import pandas as pd
from scipy.stats import norm

data = pd.DataFrame([
    {'is_call': True, 'S': 110, 'K': 100, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': False, 'S': 110, 'K': 100, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': True, 'S': 100, 'K': 100, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': False, 'S': 100, 'K': 100, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': True, 'S': 100, 'K': 110, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': False, 'S': 100, 'K': 110, 'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
])

def price(df):
    d1 = (np.log(df['S'] / df['K']) + df['T'] * (df['b'] + df['v'] ** 2 / 2)) / (df['v'] * np.sqrt(df['T']))
    d2 = d1 - df['v'] * np.sqrt(df['T'])

    return pd.Series(np.where(
        df['is_call'],
        df['S'] * np.exp((df['b'] - df['r']) * df['T']) * norm.cdf(d1) - df['K'] * np.exp(-df['r'] * df['T']) * norm.cdf(d2),
        df['K'] * np.exp(-df['r'] * df['T']) * norm.cdf(-d2) - df['S'] * np.exp((df['b'] - df['r']) * df['T']) * norm.cdf(-d1)
    ), name='price')

data['b'] = data['r'] - data['q']
x = price(data)
```

## Contributions

Contributions are welcome!

The goals of the project are centred on clarity and accuracy.

Valuable contributions include:

* More tests. Note that tests are code too, and have to be maintained. We should
    aim for the smallest complete set possible. Requests to delete are as
    relevant as those to add.
* More implementations. Pricing models should come with tests, and a numerical
    bumping framework.

The code is formatted with autopep8, and linted with pylint and mypy. Typing is
used throughout to help the reader.

Optimisations or integrations with other packages should be included as
examples only.
