import React from 'react'

import FormControlLabel from '@mui/material/FormControlLabel'
import FormControl from '@mui/material/FormControl'
import FormLabel from '@mui/material/FormLabel'
import Radio from '@mui/material/Radio'
import RadioGroup from '@mui/material/RadioGroup'

import { FieldProps } from './types'

export interface RadioFieldProps extends FieldProps {
  onChange: (value: string) => void
  options: { label: string; value: string }[]
  value: string
  row: boolean
}

const RadioField: React.FC<RadioFieldProps> = ({
  onChange,
  options,
  value,
  label,
  row
}) => (
  <FormControl>
    <FormLabel>{label}</FormLabel>
    <RadioGroup
      row={row}
      value={value}
      onChange={event => onChange(event.target.value)}
    >
      {options.map(({ label, value }) => (
        <FormControlLabel
          key={value}
          value={value}
          control={<Radio />}
          label={label}
        />
      ))}
    </RadioGroup>
  </FormControl>
)

export default RadioField
