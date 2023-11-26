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

import {
  FieldDefinition,
  NumberFieldDefinition,
  BooleanFieldDefinition
} from '../types'

export interface OptionRunnerProps {
  fields: FieldDefinition[]
  priceArgs: string[]
  analyticsArgs: Record<string, string[] | null>
  analyticImportPath: string
  numericImportPath: string
}

const OptionRunner: React.FC<OptionRunnerProps> = ({
  fields,
  priceArgs,
  analyticsArgs,
  analyticImportPath,
  numericImportPath
}) => {
  const [args, setArgs] = useState<
    Record<string, number | boolean | undefined>
  >({})
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

  const fieldProps: FieldProps[] = fields.map(fieldDefinition =>
    fieldDefinition.type === 'number'
      ? toNumberFieldProps(fieldDefinition as NumberFieldDefinition)
      : toRadioSwitchProps(fieldDefinition as BooleanFieldDefinition)
  )

  useEffect(() => {
    setArgs(
      fields.reduce(
        (obj, field) => ({ ...obj, [field.field]: field.defaultValue }),
        {}
      )
    )
  }, [fields])

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
  }, [
    pyodide,
    isRequirementsLoaded,
    args,
    analyticImportPath,
    numericImportPath,
    analyticsArgs,
    priceArgs
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

      <DynamicForm fields={fieldProps} direction="row" />

      {greeks && <OptionResultView optionResults={greeks} />}
    </Stack>
  )
}

export default OptionRunner
