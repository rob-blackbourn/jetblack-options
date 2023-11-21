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
import { runBlackScholes73 } from './blackScholes73Runner'

import type { OptionResults } from './types'

export interface BlackScholes73Props {
  pyodide: PyodideInterface
}

const toOptionalNumber = (value: string | undefined) =>
  value ? Number.parseFloat(value) : undefined

const BlackScholes73: React.FC<BlackScholes73Props> = ({ pyodide }) => {
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

    try {
      const optionResults = runBlackScholes73(
        pyodide,
        optionType === 'call',
        assetPrice,
        strikePrice,
        timeToExpiry,
        riskFreeRate,
        volatility
      )
      setGreeks(optionResults)
    } catch (error) {
      console.log(error)
    }
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

export default BlackScholes73
