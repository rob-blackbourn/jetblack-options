import { NumberFieldDefinition, BooleanFieldDefinition } from '../types'

import { ModelRoute } from './types'

export const coxRossRubenstein: ModelRoute = {
  path: '/cox-ross-rubenstein',
  name: 'Cox-Ross-Rubenstein',
  description: 'Cox-Ross-Rubenstein binomial tree',
  fields: [
    {
      label: 'Option Style',
      field: 'isEuropean',
      type: 'boolean',
      trueOption: 'European',
      falseOption: 'American',
      defaultValue: true
    } as BooleanFieldDefinition,
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
      label: 'Cost of Carry',
      field: 'costOfCarry',
      type: 'number',
      defaultValue: 0.002
    } as NumberFieldDefinition,
    {
      label: 'Volatility',
      field: 'volatility',
      type: 'number',
      defaultValue: 0.25
    } as NumberFieldDefinition,
    {
      label: 'Steps',
      field: 'steps',
      type: 'number',
      defaultValue: 200
    } as NumberFieldDefinition
  ],
  analyticImportPath: 'jetblack_options.trees.cox_ross_rubinstein',
  pricePrototype: [
    'isEuropean',
    'isCall',
    'assetPrice',
    'strikePrice',
    'timeToExpiry',
    'riskFreeRate',
    'costOfCarry',
    'volatility',
    'steps'
  ],
  greeksPrototypes: {
    delta: null,
    gamma: null,
    theta: null,
    vega: null,
    rho: null
  },
  bumpFactoryPrototype: ['isEuropean', 'isCall', 'steps'],
  bumpPrototype: [
    'assetPrice',
    'strikePrice',
    'timeToExpiry',
    'riskFreeRate',
    'costOfCarry',
    'volatility'
  ]
}

export default coxRossRubenstein
