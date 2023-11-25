import React from 'react'

import { FieldProps } from './types'

import AutoRadioField, { RadioFieldProps } from './RadioField'
import AutoTextField, { NumberFieldProps } from './NumberField'

const AutoField: React.FC<FieldProps> = props =>
  props.type === 'number' ? (
    <AutoTextField {...(props as NumberFieldProps)} />
  ) : (
    <AutoRadioField {...(props as RadioFieldProps)} />
  )

export default AutoField
