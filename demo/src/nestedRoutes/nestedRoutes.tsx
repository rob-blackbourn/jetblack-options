import { Route } from 'react-router-dom'

import Home from '../pages/Home'
import NotFound from '../pages/NotFound'

import OptionRunner from '../components/OptionRunner'

import { ModelRoute } from './types'

import blackScholes73 from './blackScholes73'
import black76 from './black76'
import blackScholesMerton from './blackScholesMerton'
import garmanKohlhagen from './garmanKohlhagen'
import generalisedBlackScholes from './generalisedBlackScholes'
import baroneAdesiWhaley from './baroneAdesiWhaley'

export const routes: ModelRoute[] = [
  blackScholes73,
  black76,
  blackScholesMerton,
  garmanKohlhagen,
  generalisedBlackScholes,
  baroneAdesiWhaley
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
        numericImportPath,
        description
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
              description={description}
            />
          }
        />
      )
    )}
    <Route path="*" element={<NotFound />} />
  </>
)
