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
  pricePrototype: string[]
  analyticsArgs: Record<string, string[] | null>
  analyticImportPath: string
  bumpFactoryPrototype: string[]
  bumpPrototype: string[]
  name: string
  description: string
}

const OptionRunner: React.FC<OptionRunnerProps> = ({
  fields,
  pricePrototype,
  analyticsArgs,
  analyticImportPath,
  bumpFactoryPrototype,
  bumpPrototype,
  name,
  description
}) => {
  const [args, setArgs] = useState<
    Record<string, number | boolean | undefined>
  >(
    fields.reduce(
      (obj, field) => ({ ...obj, [field.field]: field.defaultValue }),
      {}
    )
  )
  const [greeks, setGreeks] = useState<OptionResults>()
  const [fieldProps, setFieldProps] = useState<FieldProps[]>([])
  const { pyodide, isRequirementsLoaded } = useContext(PyodideContext)

  console.log({ fields, args })

  useEffect(() => {
    if (!args) {
      return
    }

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
    setFieldProps(fieldProps)
  }, [fields, args])

  useEffect(() => {
    if (!args) {
      return
    }

    if (
      !(
        pyodide &&
        isRequirementsLoaded &&
        Object.values(args).every(x => x != null)
      )
    ) {
      console.log('Clearing greeks')
      setGreeks(undefined)
      return
    }

    try {
      console.log('Setting greeks')
      const optionResults = valueOption(
        pyodide,
        args,
        analyticImportPath,
        pricePrototype,
        analyticsArgs,
        bumpFactoryPrototype,
        bumpPrototype
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
    analyticsArgs,
    bumpFactoryPrototype,
    bumpPrototype,
    pricePrototype
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
        <Typography variant="h4">{name}</Typography>
        <Typography variant="body1">{description}</Typography>
      </Box>

      <DynamicForm fields={fieldProps} direction="row" />

      {greeks && <OptionResultView optionResults={greeks} />}
    </Stack>
  )
}

export default OptionRunner
