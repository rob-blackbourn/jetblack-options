import React, { useContext, useEffect, useState } from 'react'

import Box from '@mui/material/Box'
import CircularProgress from '@mui/material/CircularProgress'
import Stack from '@mui/material/Stack'
import Typography from '@mui/material/Typography'

import DynamicForm, {
  FieldProps,
  NumberFieldProps,
  RadioFieldProps
} from './DynamicForm'

import OptionResultView from './OptionResultView'
import { valueOption } from './optionValuer'

import { PyodideContext } from './PythonContext'
import type { OptionResults } from './types'

export interface Black76Props {}

const Black76: React.FC<Black76Props> = () => {
  const [assetPrice, setAssetPrice] = useState<number | undefined>(100)
  const [strikePrice, setStrikePrice] = useState<number | undefined>(100)
  const [timeToExpiry, setTimeToExpiry] = useState<number | undefined>(0.5)
  const [riskFreeRate, setRiskFreeRate] = useState<number | undefined>(0.05)
  const [volatility, setVolatility] = useState<number | undefined>(0.25)
  const [optionType, setOptionType] = React.useState<'call' | 'put'>('call')
  const [greeks, setGreeks] = useState<OptionResults>()
  const { pyodide, isRequirementsLoaded } = useContext(PyodideContext)

  const handleSetOptionType = (value: string) => {
    setOptionType(value as 'call' | 'put')
  }

  const fields: FieldProps[] = [
    {
      label: 'Option Type',
      type: 'radio',
      onChange: handleSetOptionType,
      options: [
        { label: 'Call', value: 'call' },
        { label: 'Put', value: 'put' }
      ],
      value: optionType,
      row: true
    } as RadioFieldProps,
    {
      label: 'Asset Price',
      type: 'number',
      value: assetPrice,
      onChange: setAssetPrice,
      width: 200
    } as NumberFieldProps,
    {
      label: 'Strike Price',
      type: 'number',
      value: strikePrice,
      onChange: setStrikePrice,
      width: 200
    } as NumberFieldProps,
    {
      label: 'Time To Expiry',
      type: 'number',
      value: timeToExpiry,
      onChange: setTimeToExpiry,
      width: 200
    } as NumberFieldProps,
    {
      label: 'Risk Free Rate',
      type: 'number',
      value: riskFreeRate,
      onChange: setRiskFreeRate,
      width: 200
    } as NumberFieldProps,
    {
      label: 'Volatility',
      type: 'number',
      value: volatility,
      onChange: setVolatility,
      width: 200
    } as NumberFieldProps
  ]

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
        volatility
      )
    ) {
      setGreeks(undefined)
      return
    }

    try {
      const optionResults = valueOption(
        pyodide,
        {
          is_call: optionType === 'call',
          S: assetPrice,
          K: strikePrice,
          T: timeToExpiry,
          r: riskFreeRate,
          v: volatility
        },
        'jetblack_options.european.black_76',
        'jetblack_options.numeric_greeks.without_carry',
        ['is_call', 'S', 'K', 'T', 'r', 'v'],
        {
          delta: ['is_call', 'S', 'K', 'T', 'r', 'v'],
          gamma: ['S', 'K', 'T', 'r', 'v'],
          theta: ['is_call', 'S', 'K', 'T', 'r', 'v'],
          vega: ['S', 'K', 'T', 'r', 'v'],
          rho: ['is_call', 'S', 'K', 'T', 'r', 'v']
        }
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
      <Box>
        <Typography variant="body1">
          The Black model for European options on futures, bonds, and rates.
        </Typography>
      </Box>

      <DynamicForm fields={fields} direction="row" />

      {greeks && <OptionResultView optionResults={greeks} />}
    </Stack>
  )
}

export default Black76
