import componentTypes from '@data-driven-forms/react-form-renderer/component-types'
import validatorTypes from '@data-driven-forms/react-form-renderer/validator-types'

import { ModelRoute } from './types'

export const generalisedBlackScholes: ModelRoute = {
  path: '/generalised-black-scholes',
  schema: {
    title: 'Generalised Black-Scholes',
    description: 'The generalised Black-Scholes model.',
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
  analyticsPrototypes: {
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
