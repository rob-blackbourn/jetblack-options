import { PyodideInterface } from 'pyodide'

import { OptionResults } from './types'

export function runGarmanKolhagen(
  pyodide: PyodideInterface,
  is_call: boolean,
  S: number,
  K: number,
  T: number,
  r: number,
  rf: number,
  v: number
): OptionResults {
  const locals = pyodide.toPy({
    args: {
      is_call,
      S,
      K,
      T,
      r,
      rf,
      v
    }
  })

  const results = pyodide.runPython(
    `
from jetblack_options.european.garman_kolhagen import price
from jetblack_options.numeric_greeks.with_dividend_yield import NumericGreeks

is_call = args['is_call']
S = args['S']
K = args['K']
T = args['T']
r = args['r']
rf = args['rf']
v = args['v']

ng = NumericGreeks(price)
{
    'price': price(is_call, S, K, T, r, rf, v),
    'analytic': {
        'delta': None,
        'gamma': None,
        'theta': None,
        'vega': None,
        'rho': None,
    },
    'numeric': {
        'delta': ng.delta(is_call, S, K, T, r, rf, v),
        'gamma': ng.gamma(is_call, S, K, T, r, rf, v),
        'theta': ng.theta(is_call, S, K, T, r, rf, v),
        'vega': ng.vega(is_call, S, K, T, r, rf, v),
        'rho': ng.rho(is_call, S, K, T, r, rf, v),
    },
}
    `,
    { locals }
  )

  const dct = results.toJs()
  const analytic = dct.get('analytic')
  const numeric = dct.get('numeric')
  const optionResults: OptionResults = {
    price: dct.get('price'),
    analytic: {
      delta: analytic.get('delta'),
      gamma: analytic.get('gamma'),
      theta: analytic.get('theta'),
      vega: analytic.get('vega'),
      rho: analytic.get('rho')
    },
    numeric: {
      delta: numeric.get('delta'),
      gamma: numeric.get('gamma'),
      theta: numeric.get('theta'),
      vega: numeric.get('vega'),
      rho: numeric.get('rho')
    }
  }

  return optionResults
}
