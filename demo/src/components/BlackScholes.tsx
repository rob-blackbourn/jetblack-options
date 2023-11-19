import React, { useState, useEffect } from 'react'

import { PyodideInterface } from 'pyodide'

import Radio from '@mui/material/Radio'
import RadioGroup from '@mui/material/RadioGroup'
import FormControlLabel from '@mui/material/FormControlLabel'
import FormControl from '@mui/material/FormControl'
import FormLabel from '@mui/material/FormLabel'
import Stack from '@mui/material/Stack'
import TextField from '@mui/material/TextField'

interface Greeks {
  price: number
  delta: number
  gamma: number
  theta: number
  vega: number
  rho: number
}

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
  const [greeks, setGreeks] = useState<Greeks>()

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

is_call = args['is_call']
S = args['S']
K = args['K']
T = args['T']
r = args['r']
v = args['v']

{
    'price': price(is_call, S, K, T, r, v),
    'delta': delta(is_call, S, K, T, r, v),
    'gamma': gamma(S, K, T, r, v),
    'theta': theta(is_call, S, K, T, r, v),
    'vega': vega(S, K, T, r, v),
    'rho': rho(is_call, S, K, T, r, v),
}
    `,
        { locals }
      )
      .then(results => {
        const dct = results.toJs()
        setGreeks({
          price: dct.get('price'),
          delta: dct.get('delta'),
          gamma: dct.get('gamma'),
          theta: dct.get('theta'),
          vega: dct.get('vega'),
          rho: dct.get('rho')
        })
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
        {greeks && (
          <>
            <TextField
              label="Price"
              type="number"
              value={greeks.price}
              InputProps={{
                readOnly: true
              }}
              sx={{ width: 200 }}
            />
            <TextField
              label="Delta"
              type="number"
              value={greeks.delta}
              InputProps={{
                readOnly: true
              }}
              sx={{ width: 200 }}
            />
            <TextField
              label="Gamma"
              type="number"
              value={greeks.gamma}
              InputProps={{
                readOnly: true
              }}
              sx={{ width: 200 }}
            />
            <TextField
              label="Theta"
              type="number"
              value={greeks.theta}
              InputProps={{
                readOnly: true
              }}
              sx={{ width: 200 }}
            />
            <TextField
              label="Vega"
              type="number"
              value={greeks.vega}
              InputProps={{
                readOnly: true
              }}
              sx={{ width: 200 }}
            />
            <TextField
              label="Rho"
              type="number"
              value={greeks.rho}
              InputProps={{
                readOnly: true
              }}
              sx={{ width: 200 }}
            />
          </>
        )}
      </Stack>
    </Stack>
  )
}

export default BlackScholes
