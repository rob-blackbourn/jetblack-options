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

interface FieldDefinition {
  label: string
  field: string
  type: 'number' | 'boolean'
  defaultValue: boolean | number | undefined
}

interface BooleanFieldDefinition extends FieldDefinition {
  trueOption: string
  falseOption: string
  defaultValue: boolean
}

interface NumberFieldDefinition extends FieldDefinition {
  defaultValue: number
}

const fieldDefinitions: FieldDefinition[] = [
  {
    label: 'Option Type',
    field: 'isCall',
    type: 'boolean',
    trueOption: 'Call',
    falseOption: 'Put',
    defaultValue: true
  } as BooleanFieldDefinition,
  {
    label: 'Asset Price',
    field: 'assetPrice',
    type: 'number',
    defaultValue: 100
  } as NumberFieldDefinition,
  {
    label: 'Strike Price',
    field: 'strikePrice',
    type: 'number',
    defaultValue: 100
  } as NumberFieldDefinition,
  {
    label: 'Time To Expiry',
    field: 'timeToExpiry',
    type: 'number',
    defaultValue: 0.5
  } as NumberFieldDefinition,
  {
    label: 'Risk Free Rate',
    field: 'riskFreeRate',
    type: 'number',
    defaultValue: 0.005
  } as NumberFieldDefinition,
  {
    label: 'Volatility',
    field: 'volatility',
    type: 'number',
    defaultValue: 0.25
  } as NumberFieldDefinition
]

const priceArgs = [
  'isCall',
  'assetPrice',
  'strikePrice',
  'timeToExpiry',
  'riskFreeRate',
  'volatility'
]

const analyticsArgs = {
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

const analyticImportPath = 'jetblack_options.european.black_76'
const numericImportPath = 'jetblack_options.numeric_greeks.without_carry'

const Black76: React.FC<Black76Props> = () => {
  const [args, setArgs] = useState<
    Record<string, number | boolean | undefined>
  >(
    fieldDefinitions.reduce(
      (obj, field) => ({ ...obj, [field.field]: field.defaultValue }),
      {}
    )
  )
  const [greeks, setGreeks] = useState<OptionResults>()
  const { pyodide, isRequirementsLoaded } = useContext(PyodideContext)

  const toNumberFieldProps = ({
    label,
    field
  }: NumberFieldDefinition): NumberFieldProps => ({
    label,
    type: 'number',
    onChange: (value: number | undefined) =>
      setArgs(state => ({ ...state, [field]: value })),
    value: args[field] as number | undefined,
    width: 200
  })

  const toRadioSwitchProps = ({
    label,
    field,
    trueOption,
    falseOption
  }: BooleanFieldDefinition): RadioSwitchFieldProps => ({
    label,
    type: 'radio-switch',
    trueOption,
    falseOption,
    onChange: (value: boolean | undefined) =>
      setArgs(state => ({ ...state, [field]: value })),
    value: args[field] as boolean,
    row: true
  })

  const fieldProps: FieldProps[] = fieldDefinitions.map(fieldDefinition =>
    fieldDefinition.type === 'number'
      ? toNumberFieldProps(fieldDefinition as NumberFieldDefinition)
      : toRadioSwitchProps(fieldDefinition as BooleanFieldDefinition)
  )

  useEffect(() => {
    if (
      !(
        pyodide &&
        isRequirementsLoaded &&
        args.assetPrice &&
        args.strikePrice &&
        args.timeToExpiry &&
        args.riskFreeRate &&
        args.volatility &&
        !Object.values(args).some(x => x === undefined)
      )
    ) {
      setGreeks(undefined)
      return
    }

    try {
      const optionResults = valueOption(
        pyodide,
        args,
        analyticImportPath,
        numericImportPath,
        priceArgs,
        analyticsArgs
      )
      setGreeks(optionResults)
    } catch (error) {
      console.log(error)
    }
  }, [pyodide, isRequirementsLoaded, args])

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

      <DynamicForm fields={fieldProps} direction="row" />

      {greeks && <OptionResultView optionResults={greeks} />}
    </Stack>
  )
}

export default Black76
