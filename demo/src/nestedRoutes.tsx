import { Route } from 'react-router-dom'

import Home from './pages/Home'
import NotFound from './pages/NotFound'

import OptionRunner from './components/OptionRunner'

import {
  FieldDefinition,
  NumberFieldDefinition,
  BooleanFieldDefinition
} from './types'

export interface ModelRoute {
  path: string
  name: string
  fields: FieldDefinition[]
  analyticImportPath: string
  numericImportPath: string
  pricePrototype: string[]
  greeksPrototypes: Record<string, string[] | null>
}

export const routes: ModelRoute[] = [
  {
    path: '/black-scholes-73',
    name: 'Black-Scholes 73',
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
    numericImportPath: 'jetblack_options.numeric_greeks.without_carry',
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
    }
  },
  {
    path: '/black-76',
    name: 'Black 76',
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
    analyticImportPath: 'jetblack_options.european.black_76',
    numericImportPath: 'jetblack_options.numeric_greeks.without_carry',
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
    }
  },
  {
    path: '/black-scholes-merton',
    name: 'Black-Scholes-Merton',
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
        label: 'Dividend Yield',
        field: 'dividendYield',
        type: 'number',
        defaultValue: 0.002
      } as NumberFieldDefinition,
      {
        label: 'Volatility',
        field: 'volatility',
        type: 'number',
        defaultValue: 0.25
      } as NumberFieldDefinition
    ],
    analyticImportPath: 'jetblack_options.european.black_scholes_merton',
    numericImportPath: 'jetblack_options.numeric_greeks.with_dividend_yield',
    pricePrototype: [
      'isCall',
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'dividendYield',
      'volatility'
    ],
    greeksPrototypes: {
      delta: [
        'isCall',
        'assetPrice',
        'strikePrice',
        'timeToExpiry',
        'riskFreeRate',
        'dividendYield',
        'volatility'
      ],
      gamma: [
        'assetPrice',
        'strikePrice',
        'timeToExpiry',
        'riskFreeRate',
        'dividendYield',
        'volatility'
      ],
      theta: [
        'isCall',
        'assetPrice',
        'strikePrice',
        'timeToExpiry',
        'riskFreeRate',
        'dividendYield',
        'volatility'
      ],
      vega: [
        'assetPrice',
        'strikePrice',
        'timeToExpiry',
        'riskFreeRate',
        'dividendYield',
        'volatility'
      ],
      rho: [
        'isCall',
        'assetPrice',
        'strikePrice',
        'timeToExpiry',
        'riskFreeRate',
        'dividendYield',
        'volatility'
      ]
    }
  },
  {
    path: '/garman-kohlhagen',
    name: 'Garman Kohlhagen',
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
        label: 'Base Yield',
        field: 'riskFreeRate',
        type: 'number',
        defaultValue: 0.005
      } as NumberFieldDefinition,
      {
        label: 'Quote Yield',
        field: 'quoteRiskFreeRate',
        type: 'number',
        defaultValue: 0.002
      } as NumberFieldDefinition,
      {
        label: 'Volatility',
        field: 'volatility',
        type: 'number',
        defaultValue: 0.25
      } as NumberFieldDefinition
    ],
    analyticImportPath: 'jetblack_options.european.garman_kolhagen',
    numericImportPath: 'jetblack_options.numeric_greeks.with_dividend_yield',
    pricePrototype: [
      'isCall',
      'assetPrice',
      'strikePrice',
      'timeToExpiry',
      'riskFreeRate',
      'quoteRiskFreeRate',
      'volatility'
    ],
    greeksPrototypes: {
      delta: null,
      gamma: null,
      theta: null,
      vega: null,
      rho: null
    }
  },
  {
    path: '/generalised-black-scholes',
    name: 'Generalised Black-Scholes',
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
        label: 'Carry RateYield',
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
    }
  },
  {
    path: '/baron-adesi-whaley',
    name: 'Baron-Adesi-Whaley',
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
        label: 'Carry RateYield',
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
    analyticImportPath: 'jetblack_options.american.barone_adesi_whaley',
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
]

export const renderNestedRoutes = () => (
  <>
    <Route index element={<Home />} />
    {routes.map(
      ({
        path,
        fields,
        pricePrototype,
        greeksPrototypes,
        analyticImportPath,
        numericImportPath
      }) => (
        <Route
          key={path}
          path={path}
          element={
            <OptionRunner
              fields={fields}
              priceArgs={pricePrototype}
              analyticsArgs={greeksPrototypes}
              analyticImportPath={analyticImportPath}
              numericImportPath={numericImportPath}
            />
          }
        />
      )
    )}
    <Route path="*" element={<NotFound />} />
  </>
)
