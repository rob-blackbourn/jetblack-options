import { NumberFieldDefinition, BooleanFieldDefinition } from '../types'

import { ModelRoute } from './types'

export const bjerksundStensland2002: ModelRoute = {
  path: '/bjerksund-stensland-2002',
  name: 'Bjerksund-Stensland 2002',
  description: 'Bjerksund-Stensland 2002',
  fields: [
    {
      label: 'Option Type',
      field: 'isCall',
      type: 'boolean',
      trueOption: 'Call',
      falseOption: 'Put',
      defaultValue: true
    } as BooleanFieldDefinition,
    {
      label: 'Asset Price',
      field: 'assetPrice',
      type: 'number',
      defaultValue: 100
    } as NumberFieldDefinition,
    {
      label: 'Strike Price',
      field: 'strikePrice',
      type: 'number',
      defaultValue: 100
    } as NumberFieldDefinition,
    {
      label: 'Time To Expiry',
      field: 'timeToExpiry',
      type: 'number',
      defaultValue: 0.5
    } as NumberFieldDefinition,
    {
      label: 'Risk Free Rate',
      field: 'riskFreeRate',
      type: 'number',
      defaultValue: 0.005
    } as NumberFieldDefinition,
    {
      label: 'Carry Rate',
      field: 'costOfCarry',
      type: 'number',
      defaultValue: 0.003
    } as NumberFieldDefinition,
    {
      label: 'Volatility',
      field: 'volatility',
      type: 'number',
      defaultValue: 0.25
    } as NumberFieldDefinition
  ],
  analyticImportPath: 'jetblack_options.american.bjerksund_stensland_2002',
  numericImportPath: 'jetblack_options.numeric_greeks.with_carry',
  pricePrototype: [
    'isCall',
    'assetPrice',
    'strikePrice',
    'timeToExpiry',
    'riskFreeRate',
    'costOfCarry',
    'volatility'
  ],
  greeksPrototypes: {
    delta: null,
    gamma: null,
    theta: null,
    vega: null,
    rho: null
  }
}

export default bjerksundStensland2002
