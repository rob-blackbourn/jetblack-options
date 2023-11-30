import { Suspense } from 'react'

import { HashRouter as Router, Route, Routes } from 'react-router-dom'

import Loading from './components/Loading'

import Layout from './pages/Layout'
import Home from './pages/Home'
import NotFound from './pages/NotFound'
import OptionRunner from './components/OptionRunner'

import { routes } from './nestedRoutes'

function RoutedApp() {
  return (
    <Router>
      <Suspense fallback={<Loading />}>
        <Routes>
          <Route path="/" element={<Layout />}>
            <Route index element={<Home />} />
            {routes.map(
              ({
                path,
                fields,
                pricePrototype,
                analyticsPrototypes: greeksPrototypes,
                analyticImportPath,
                bumpFactoryPrototype,
                bumpPrototype,
                description
              }) => (
                <Route
                  key={path}
                  path={path}
                  element={
                    <OptionRunner
                      fields={fields}
                      pricePrototype={pricePrototype}
                      analyticsArgs={greeksPrototypes}
                      analyticImportPath={analyticImportPath}
                      bumpFactoryPrototype={bumpFactoryPrototype}
                      bumpPrototype={bumpPrototype}
                      description={description}
                    />
                  }
                />
              )
            )}
            <Route path="*" element={<NotFound />} />
          </Route>
        </Routes>
      </Suspense>
    </Router>
  )
}

export default RoutedApp
