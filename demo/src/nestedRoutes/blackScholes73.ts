import { NumberFieldDefinition, BooleanFieldDefinition } from '../types'

import { ModelRoute } from './types'

const blackScholes73: ModelRoute = {
  path: '/black-scholes-73',
  name: 'Black-Scholes 73',
  description:
    'The original Black-Scholes for European options on non-dividend paying stock.',
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
      label: 'Volatility',
      field: 'volatility',
      type: 'number',
      defaultValue: 0.25
    } as NumberFieldDefinition
  ],
  analyticImportPath: 'jetblack_options.european.black_scholes_73',
  pricePrototype: [
    'isCall',
    'assetPrice',
    'strikePrice',
    'timeToExpiry',
    'riskFreeRate',
    'volatility'
  ],
  greeksPrototypes: {
    delta: [
      'isCall',
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'volatility'
    ],
    gamma: [
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'volatility'
    ],
    theta: [
      'isCall',
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'volatility'
    ],
    vega: [
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'volatility'
    ],
    rho: [
      'isCall',
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'volatility'
    ]
  },
  bumpFactoryPrototype: ['isCall'],
  bumpPrototype: [
    'assetPrice',
    'strikePrice',
    'timeToExpiry',
    'riskFreeRate',
    'volatility'
  ]
}

export default blackScholes73
