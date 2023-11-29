import { NumberFieldDefinition, BooleanFieldDefinition } from '../types'

import { ModelRoute } from './types'

export const generalisedBlackScholes: ModelRoute = {
  path: '/generalised-black-scholes',
  name: 'Generalised Black-Scholes',
  description: 'The generalised Black-Scholes model.',
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
  analyticImportPath: 'jetblack_options.european.generalised_black_scholes',
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
    delta: [
      'isCall',
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'costOfCarry',
      'volatility'
    ],
    gamma: [
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'costOfCarry',
      'volatility'
    ],
    theta: [
      'isCall',
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'costOfCarry',
      'volatility'
    ],
    vega: [
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'costOfCarry',
      'volatility'
    ],
    rho: [
      'isCall',
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'costOfCarry',
      'volatility'
    ]
  },
  bumpFactoryPrototype: ['isCall'],
  bumpPrototype: [
    'assetPrice',
    'strikePrice',
    'timeToExpiry',
    'riskFreeRate',
    'costOfCarry',
    'volatility'
  ]
}

export default generalisedBlackScholes
