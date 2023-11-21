import React, { useContext, useEffect, useState } from 'react'

import CircularProgress from '@mui/material/CircularProgress'
import FormControlLabel from '@mui/material/FormControlLabel'
import FormControl from '@mui/material/FormControl'
import FormLabel from '@mui/material/FormLabel'
import Radio from '@mui/material/Radio'
import RadioGroup from '@mui/material/RadioGroup'
import Stack from '@mui/material/Stack'
import TextField from '@mui/material/TextField'

import OptionResultView from './OptionResultView'
import { runBlackScholesMerton } from './blackScholesMertonRunner'

import { PyodideContext } from './PythonContext'
import type { OptionResults } from './types'

export interface BlackScholesMertonProps {}

const toOptionalNumber = (value: string | undefined) =>
  value ? Number.parseFloat(value) : undefined

const BlackScholesMerton: React.FC<BlackScholesMertonProps> = () => {
  const [assetPrice, setAssetPrice] = useState<number | undefined>(100)
  const [strikePrice, setStrikePrice] = useState<number | undefined>(100)
  const [timeToExpiry, setTimeToExpiry] = useState<number | undefined>(0.5)
  const [riskFreeRate, setRiskFreeRate] = useState<number | undefined>(0.05)
  const [dividendYield, setDividendYield] = useState<number | undefined>(0.02)
  const [volatility, setVolatility] = useState<number | undefined>(0.25)
  const [optionType, setOptionType] = React.useState<'call' | 'put'>('call')
  const [greeks, setGreeks] = useState<OptionResults>()
  const { pyodide, isRequirementsLoaded } = useContext(PyodideContext)

  const handleSetOptionType = (event: React.ChangeEvent<HTMLInputElement>) => {
    setOptionType((event.target as HTMLInputElement).value as 'call' | 'put')
  }

  useEffect(() => {
    if (
      !(
        pyodide &&
        isRequirementsLoaded &&
        optionType &&
        assetPrice &&
        strikePrice &&
        timeToExpiry &&
        riskFreeRate &&
        dividendYield &&
        volatility
      )
    ) {
      setGreeks(undefined)
      return
    }

    try {
      const optionResults = runBlackScholesMerton(
        pyodide,
        optionType === 'call',
        assetPrice,
        strikePrice,
        timeToExpiry,
        riskFreeRate,
        dividendYield,
        volatility
      )
      setGreeks(optionResults)
    } catch (error) {
      console.log(error)
    }
  }, [
    pyodide,
    isRequirementsLoaded,
    optionType,
    assetPrice,
    strikePrice,
    timeToExpiry,
    riskFreeRate,
    dividendYield,
    volatility
  ])

  if (!(pyodide && isRequirementsLoaded)) {
    return (
      <Stack>
        <CircularProgress />
      </Stack>
    )
  }

  return (
    <Stack direction="column" spacing={2}>
      <Stack
        direction="row"
        spacing={1}
        alignItems="end"
        justifyContent="flex-start"
      >
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
          sx={{ width: 150 }}
        />
        <TextField
          label="Strike Price"
          type="number"
          value={strikePrice}
          onChange={event =>
            setStrikePrice(toOptionalNumber(event.target.value))
          }
          sx={{ width: 150 }}
        />
        <TextField
          label="Time to Expiry"
          type="number"
          value={timeToExpiry}
          onChange={event =>
            setTimeToExpiry(toOptionalNumber(event.target.value))
          }
          sx={{ width: 150 }}
        />
        <TextField
          label="Risk Free Rate"
          type="number"
          value={riskFreeRate}
          onChange={event =>
            setRiskFreeRate(toOptionalNumber(event.target.value))
          }
          sx={{ width: 150 }}
        />
        <TextField
          label="Dividend Yield"
          type="number"
          value={dividendYield}
          onChange={event =>
            setDividendYield(toOptionalNumber(event.target.value))
          }
          sx={{ width: 150 }}
        />
        <TextField
          label="Volatility"
          type="number"
          value={volatility}
          onChange={event =>
            setVolatility(toOptionalNumber(event.target.value))
          }
          sx={{ width: 150 }}
        />
      </Stack>
      <Stack direction="column" spacing={2}>
        {greeks && <OptionResultView optionResults={greeks} />}
      </Stack>
    </Stack>
  )
}

export default BlackScholesMerton
