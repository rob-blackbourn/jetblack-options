import React from 'react'

import TextField from '@mui/material/TextField'

import { FieldProps } from './types'

const toOptionalNumber = (value: string | undefined) =>
  value ? Number.parseFloat(value) : undefined

export interface NumberFieldProps extends FieldProps {
  label: string
  onChange: (value: number | undefined) => void
  value: number | undefined
  width: number
}

const NumberField: React.FC<NumberFieldProps> = ({
  onChange,
  label,
  value,
  width
}) => (
  <TextField
    label={label}
    type="number"
    value={value}
    onChange={event => onChange(toOptionalNumber(event.target.value))}
    sx={{ width }}
  />
)

export default NumberField
