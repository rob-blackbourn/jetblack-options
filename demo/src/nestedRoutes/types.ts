import { FieldDefinition } from '../types'

export interface ModelRoute {
  path: string
  name: string
  description: string
  fields: FieldDefinition[]
  analyticImportPath: string
  numericImportPath: string
  pricePrototype: string[]
  greeksPrototypes: Record<string, string[] | null>
}
