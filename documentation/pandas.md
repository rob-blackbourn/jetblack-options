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

There are some more examples in the `scratch` folder.
