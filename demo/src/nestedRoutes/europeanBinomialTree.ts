import componentTypes from '@data-driven-forms/react-form-renderer/component-types'
import validatorTypes from '@data-driven-forms/react-form-renderer/validator-types'

import { ModelRoute } from './types'

export const europeanBinomialTree: ModelRoute = {
  path: '/european-binomial-tree',
  schema: {
    title: 'European Binomial Tree',
    description:
      'A binomial tree model for pricing European options using cost of carry.',
    fields: [
      {
        name: 'isCall',
        label: 'Option Type',
        component: componentTypes.SWITCH,
        onText: 'Call',
        offText: 'Put',
        initialValue: true,
        dataType: 'boolean'
      },
      {
        name: 'assetPrice',
        label: 'Asset Price',
        component: componentTypes.TEXT_FIELD,
        type: 'number',
        dataType: 'number',
        validate: [
          {
            type: validatorTypes.REQUIRED
          }
        ],
        initialValue: 100
      },
      {
        name: 'strikePrice',
        label: 'Strike Price',
        component: componentTypes.TEXT_FIELD,
        type: 'number',
        dataType: 'number',
        validate: [
          {
            type: validatorTypes.REQUIRED
          }
        ],
        initialValue: 100
      },
      {
        name: 'timeToExpiry',
        label: 'Time To Expiry',
        component: componentTypes.TEXT_FIELD,
        type: 'number',
        dataType: 'number',
        validate: [
          {
            type: validatorTypes.REQUIRED
          }
        ],
        initialValue: 0.5
      },
      {
        name: 'riskFreeRate',
        label: 'Risk Free Rate',
        component: componentTypes.TEXT_FIELD,
        type: 'number',
        dataType: 'number',
        validate: [
          {
            type: validatorTypes.REQUIRED
          }
        ],
        initialValue: 0.005
      },
      {
        name: 'costOfCarry',
        label: 'Time To Expiry',
        component: componentTypes.TEXT_FIELD,
        type: 'number',
        dataType: 'number',
        validate: [
          {
            type: validatorTypes.REQUIRED
          }
        ],
        initialValue: 0.003
      },
      {
        name: 'volatility',
        label: 'Volatility',
        component: componentTypes.TEXT_FIELD,
        type: 'number',
        dataType: 'number',
        validate: [
          {
            type: validatorTypes.REQUIRED
          }
        ],
        initialValue: 0.25
      }
    ]
  },
  analyticImportPath: 'jetblack_options.trees.european_binomial',
  pricePrototype: [
    'isCall',
    'assetPrice',
    'strikePrice',
    'timeToExpiry',
    'riskFreeRate',
    'costOfCarry',
    'volatility',
    'steps'
  ],
  analyticsPrototypes: {
    delta: null,
    gamma: null,
    theta: null,
    vega: null,
    rho: null
  },
  bumpFactoryPrototype: ['isCall', 'steps'],
  bumpPrototype: [
    'assetPrice',
    'strikePrice',
    'timeToExpiry',
    'riskFreeRate',
    'costOfCarry',
    'volatility'
  ]
}

export default europeanBinomialTree
