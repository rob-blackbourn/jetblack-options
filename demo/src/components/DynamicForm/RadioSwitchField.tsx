import React from 'react'

import FormControlLabel from '@mui/material/FormControlLabel'
import FormControl from '@mui/material/FormControl'
import FormLabel from '@mui/material/FormLabel'
import Radio from '@mui/material/Radio'
import RadioGroup from '@mui/material/RadioGroup'

import { FieldProps } from './types'

export interface RadioSwitchFieldProps extends FieldProps {
  onChange: (value: boolean) => void
  trueOption: string
  falseOption: string
  value: boolean
  row: boolean
}

const RadioSwitchField: React.FC<RadioSwitchFieldProps> = ({
  onChange,
  trueOption,
  falseOption,
  value,
  label,
  row
}) => (
  <FormControl>
    <FormLabel>{label}</FormLabel>
    <RadioGroup
      row={row}
      value={value ? 'true' : 'false'}
      onChange={event => onChange(event.target.value === 'true')}
    >
      <FormControlLabel value={'true'} control={<Radio />} label={trueOption} />
      <FormControlLabel
        value={'false'}
        control={<Radio />}
        label={falseOption}
      />
    </RadioGroup>
  </FormControl>
)

export default RadioSwitchField
