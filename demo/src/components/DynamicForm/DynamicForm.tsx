import React from 'react'

import Stack from '@mui/material/Stack'

import Field from './Field'
import { FieldProps } from './types'

export interface DynamicFormProps {
  fields: FieldProps[]
  direction: 'row' | 'row-reverse' | 'column' | 'column-reverse' | undefined
}

const DynamicForm: React.FC<DynamicFormProps> = ({ fields, direction }) => (
  <Stack
    direction={direction}
    spacing={1}
    alignItems="end"
    justifyContent="flex-start"
  >
    {fields.map((field, index) => (
      <Field key={index} {...field} />
    ))}
  </Stack>
)

export default DynamicForm
