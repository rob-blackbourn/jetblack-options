import { FieldDefinition } from '../types'

export interface ModelRoute {
  path: string
  name: string
  description: string
  fields: FieldDefinition[]
  analyticImportPath: string
  pricePrototype: string[]
  analyticsPrototypes: Record<string, string[] | null>
  bumpFactoryPrototype: string[]
  bumpPrototype: string[]
}
