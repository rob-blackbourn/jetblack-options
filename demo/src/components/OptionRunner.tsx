import React, { useContext, useEffect, useState } from 'react'

import CircularProgress from '@mui/material/CircularProgress'
import Stack from '@mui/material/Stack'

import Grid from '@mui/material/Grid'
import FormRenderer from '@data-driven-forms/react-form-renderer/form-renderer'
import componentTypes from '@data-driven-forms/react-form-renderer/component-types'
import FormTemplate from '@data-driven-forms/mui-component-mapper/form-template'
import Switch from '@data-driven-forms/mui-component-mapper/switch'
import TextField from '@data-driven-forms/mui-component-mapper/text-field'
import { Schema, ComponentMapper } from '@data-driven-forms/react-form-renderer'
import { FormTemplateCommonProps } from '@data-driven-forms/common/form-template'

import OptionResultView from './OptionResultView'
import { valueOption } from './optionValuer'

import { PyodideContext } from './PythonContext'
import type { OptionResults } from './types'

export interface OptionRunnerProps {
  schema: Schema
  pricePrototype: string[]
  analyticsArgs: Record<string, string[] | null>
  analyticImportPath: string
  bumpFactoryPrototype: string[]
  bumpPrototype: string[]
}

const componentMapper: ComponentMapper = {
  [componentTypes.TEXT_FIELD]: TextField,
  [componentTypes.SWITCH]: Switch
}

const FormTemplateCanReset = (props: FormTemplateCommonProps) => (
  <FormTemplate {...props} />
)

const OptionRunner: React.FC<OptionRunnerProps> = ({
  schema,
  pricePrototype,
  analyticsArgs,
  analyticImportPath,
  bumpFactoryPrototype,
  bumpPrototype
}) => {
  const [greeks, setGreeks] = useState<OptionResults>()
  const { pyodide, isRequirementsLoaded } = useContext(PyodideContext)

  const handleSubmit = (args: object) => {
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
  }

  useEffect(() => {
    setGreeks(undefined)
  }, [schema])

  if (!(pyodide && isRequirementsLoaded)) {
    return (
      <Stack>
        <CircularProgress />
      </Stack>
    )
  }

  return (
    <Grid spacing={4} container>
      <FormRenderer
        componentMapper={componentMapper}
        FormTemplate={FormTemplateCanReset}
        schema={schema}
        onSubmit={args => handleSubmit(args)}
      />

      {greeks && (
        <OptionResultView optionResults={greeks} sx={{ marginLeft: 4 }} />
      )}
    </Grid>
  )
}

export default OptionRunner
