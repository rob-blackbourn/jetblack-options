export interface FieldDefinition {
  label: string
  field: string
  type: 'number' | 'boolean'
  defaultValue: boolean | number | undefined
}

export interface BooleanFieldDefinition extends FieldDefinition {
  trueOption: string
  falseOption: string
  defaultValue: boolean
}

export interface NumberFieldDefinition extends FieldDefinition {
  defaultValue: number
}
