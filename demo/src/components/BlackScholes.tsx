import React, { useState, useEffect } from 'react'

import { PyodideInterface } from 'pyodide'

import Radio from '@mui/material/Radio'
import RadioGroup from '@mui/material/RadioGroup'
import FormControlLabel from '@mui/material/FormControlLabel'
import FormControl from '@mui/material/FormControl'
import FormLabel from '@mui/material/FormLabel'
import Stack from '@mui/material/Stack'
import TextField from '@mui/material/TextField'

import OptionResultView from './OptionResultView'

import type { OptionResults } from './types'
export interface BlackScholesProps {
  pyodide: PyodideInterface
}

const toOptionalNumber = (value: string | undefined) =>
  value ? Number.parseFloat(value) : undefined

const BlackScholes: React.FC<BlackScholesProps> = ({ pyodide }) => {
  const [assetPrice, setAssetPrice] = useState<number | undefined>(100)
  const [strikePrice, setStrikePrice] = useState<number | undefined>(100)
  const [timeToExpiry, setTimeToExpiry] = useState<number | undefined>(0.5)
  const [riskFreeRate, setRiskFreeRate] = useState<number | undefined>(0.05)
  const [volatility, setVolatility] = useState<number | undefined>(0.25)
  const [optionType, setOptionType] = React.useState<'call' | 'put'>('call')
  const [greeks, setGreeks] = useState<OptionResults>()

  const handleSetOptionType = (event: React.ChangeEvent<HTMLInputElement>) => {
    setOptionType((event.target as HTMLInputElement).value as 'call' | 'put')
  }

  useEffect(() => {
    if (
      !(
        pyodide &&
        optionType &&
        assetPrice &&
        strikePrice &&
        timeToExpiry &&
        riskFreeRate &&
        volatility
      )
    ) {
      return
    }

    const locals = pyodide.toPy({
      args: {
        is_call: optionType === 'call',
        S: assetPrice,
        K: strikePrice,
        T: timeToExpiry,
        r: riskFreeRate,
        v: volatility
      }
    })
    pyodide
      .runPythonAsync(
        `
from jetblack_options.european.black_scholes_73 import (
  price,
  delta,
  gamma,
  theta,
  vega,
  rho
)
from jetblack_options.numeric_greeks.without_carry import NumericGreeks

is_call = args['is_call']
S = args['S']
K = args['K']
T = args['T']
r = args['r']
v = args['v']

ng = NumericGreeks(price)
{
    'price': price(is_call, S, K, T, r, v),
    'analytic': {
        'delta': delta(is_call, S, K, T, r, v),
        'gamma': gamma(S, K, T, r, v),
        'theta': theta(is_call, S, K, T, r, v),
        'vega': vega(S, K, T, r, v),
        'rho': rho(is_call, S, K, T, r, v),
    },
    'numeric': {
        'delta': ng.delta(is_call, S, K, T, r, v),
        'gamma': ng.gamma(is_call, S, K, T, r, v),
        'theta': ng.theta(is_call, S, K, T, r, v),
        'vega': ng.vega(is_call, S, K, T, r, v),
        'rho': ng.rho(is_call, S, K, T, r, v),
    },
}
    `,
        { locals }
      )
      .then(results => {
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
        setGreeks(optionResults)
      })
      .catch(error => {
        console.log(error)
      })
  }, [
    pyodide,
    optionType,
    assetPrice,
    strikePrice,
    timeToExpiry,
    riskFreeRate,
    volatility
  ])

  return (
    <Stack direction="row" spacing={2}>
      <Stack direction="column" spacing={2}>
        <FormControl>
          <FormLabel>Option Type</FormLabel>
          <RadioGroup row value={optionType} onChange={handleSetOptionType}>
            <FormControlLabel value="call" control={<Radio />} label="Call" />
            <FormControlLabel value="put" control={<Radio />} label="Put" />
          </RadioGroup>
        </FormControl>
        <TextField
          label="Asset Price"
          type="number"
          value={assetPrice}
          onChange={event =>
            setAssetPrice(toOptionalNumber(event.target.value))
          }
          sx={{ width: 200 }}
        />
        <TextField
          label="Strike Price"
          type="number"
          value={strikePrice}
          onChange={event =>
            setStrikePrice(toOptionalNumber(event.target.value))
          }
          sx={{ width: 200 }}
        />
        <TextField
          label="Time to Expiry"
          type="number"
          value={timeToExpiry}
          onChange={event =>
            setTimeToExpiry(toOptionalNumber(event.target.value))
          }
          sx={{ width: 200 }}
        />
        <TextField
          label="Risk Free Rate"
          type="number"
          value={riskFreeRate}
          onChange={event =>
            setRiskFreeRate(toOptionalNumber(event.target.value))
          }
          sx={{ width: 200 }}
        />
        <TextField
          label="Volatility"
          type="number"
          value={volatility}
          onChange={event =>
            setVolatility(toOptionalNumber(event.target.value))
          }
          sx={{ width: 200 }}
        />
      </Stack>
      <Stack direction="column" spacing={2}>
        {greeks && <OptionResultView optionResults={greeks} />}
      </Stack>
    </Stack>
  )
}

export default BlackScholes
