Variants of Black scholes can be used to price: stocks with dividends, futures,
and currencies. The *generalized* Black Scholes formulation unifies these three.
variants into a single formula.

Instead of taking a dividend yield, or foreign currency rate, the formula takes
a *cost of carry* (typically represented by `b`), where:

* $b = r$ for European options on non-dividend paying stock.
* $b = r - q$ for European options and indices on dividend paying stock.
* $b = 0$ for pricing European futures options.
* $b = r - r_f$ for pricing European currency options.

The *cost of carry* methodology is not restricted to the Black Scholes, and is
used in many of the formula in this package.

## What next ?

[Using pandas](./pandas.md)