from jetblack_options.european.black_scholes_merton import (
    price,
    vega
)

v = vega(100.0, 100.0, 365.0, 0.05, 0.01, 0.25)
