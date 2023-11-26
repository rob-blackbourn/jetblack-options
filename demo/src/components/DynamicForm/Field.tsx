import React from 'react'

import { FieldProps } from './types'

import RadioField, { RadioFieldProps } from './RadioField'
import RadioSwitchField, { RadioSwitchFieldProps } from './RadioSwitchField'
import NumberField, { NumberFieldProps } from './NumberField'

const Field: React.FC<FieldProps> = props =>
  props.type === 'number' ? (
    <NumberField {...(props as NumberFieldProps)} />
  ) : props.type === 'radio' ? (
    <RadioField {...(props as RadioFieldProps)} />
  ) : (
    <RadioSwitchField {...(props as RadioSwitchFieldProps)} />
  )

export default Field
