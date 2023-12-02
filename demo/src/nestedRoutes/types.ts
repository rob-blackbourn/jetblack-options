import { Schema } from '@data-driven-forms/react-form-renderer'

export interface ModelRoute {
  path: string
  schema: Schema
  analyticImportPath: string
  pricePrototype: string[]
  analyticsPrototypes: Record<string, string[] | null>
  bumpFactoryPrototype: string[]
  bumpPrototype: string[]
}
