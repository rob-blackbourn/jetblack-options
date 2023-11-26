import React, { useContext, useEffect, useState } from 'react'

import Box from '@mui/material/Box'
import CircularProgress from '@mui/material/CircularProgress'
import Stack from '@mui/material/Stack'
import Typography from '@mui/material/Typography'

import DynamicForm, {
  FieldProps,
  NumberFieldProps,
  RadioSwitchFieldProps
} from './DynamicForm'

import OptionResultView from './OptionResultView'
import { valueOption } from './optionValuer'

import { PyodideContext } from './PythonContext'
import type { OptionResults } from './types'

export interface Black76Props {}

const Black76: React.FC<Black76Props> = () => {
  const [isCall, setIsCall] = React.useState(true)
  const [assetPrice, setAssetPrice] = useState<number | undefined>(100)
  const [strikePrice, setStrikePrice] = useState<number | undefined>(100)
  const [timeToExpiry, setTimeToExpiry] = useState<number | undefined>(0.5)
  const [riskFreeRate, setRiskFreeRate] = useState<number | undefined>(0.05)
  const [volatility, setVolatility] = useState<number | undefined>(0.25)
  const [greeks, setGreeks] = useState<OptionResults>()
  const { pyodide, isRequirementsLoaded } = useContext(PyodideContext)

  const handleSetIsCall = (value: boolean) => {
    setIsCall(value)
  }

  const fields: FieldProps[] = [
    {
      label: 'Option Type',
      type: 'radio-switch',
      onChange: handleSetIsCall,
      trueOption: 'Call',
      falseOption: 'Put',
      value: isCall,
      row: true
    } as RadioSwitchFieldProps,
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
          isCall,
          assetPrice,
          strikePrice,
          timeToExpiry,
          riskFreeRate,
          volatility
        },
        'jetblack_options.european.black_76',
        'jetblack_options.numeric_greeks.without_carry',
        [
          'isCall',
          'assetPrice',
          'strikePrice',
          'timeToExpiry',
          'riskFreeRate',
          'volatility'
        ],
        {
          delta: [
            'isCall',
            'assetPrice',
            'strikePrice',
            'timeToExpiry',
            'riskFreeRate',
            'volatility'
          ],
          gamma: [
            'assetPrice',
            'strikePrice',
            'timeToExpiry',
            'riskFreeRate',
            'volatility'
          ],
          theta: [
            'isCall',
            'assetPrice',
            'strikePrice',
            'timeToExpiry',
            'riskFreeRate',
            'volatility'
          ],
          vega: [
            'assetPrice',
            'strikePrice',
            'timeToExpiry',
            'riskFreeRate',
            'volatility'
          ],
          rho: [
            'isCall',
            'assetPrice',
            'strikePrice',
            'timeToExpiry',
            'riskFreeRate',
            'volatility'
          ]
        }
      )
      setGreeks(optionResults)
    } catch (error) {
      console.log(error)
    }
  }, [
    pyodide,
    isRequirementsLoaded,
    isCall,
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
