import { lazy, LazyExoticComponent, FC } from 'react'
import { Route } from 'react-router-dom'

import Home from './pages/Home'
import NotFound from './pages/NotFound'

export interface ModelRoute {
  path: string
  name: string
  element: LazyExoticComponent<FC<object>>
}

export const routes: ModelRoute[] = [
  {
    path: '/black-scholes-73',
    name: 'Black-Scholes 73',
    element: lazy(() => import('./components/BlackScholes73'))
  },
  {
    path: '/black-76',
    name: 'Black 76',
    element: lazy(() => import('./components/Black76'))
  },
  {
    path: '/black-scholes-merton',
    name: 'Black-Scholes-Merton',
    element: lazy(() => import('./components/BlackScholes73'))
  },
  {
    path: '/garman-kohlhagen',
    name: 'Garman Kohlhagen',
    element: lazy(() => import('./components/GarmanKohlhagen'))
  },
  {
    path: '/generalised-black-scholes',
    name: 'Generalised Black-Scholes',
    element: lazy(() => import('./components/GeneralisedBlackScholes'))
  },
  {
    path: '/baron-adesi-whaley',
    name: 'Baron-Adesi-Whaley',
    element: lazy(() => import('./components/BaroneAdesiWhaley'))
  }
]

export const renderNestedRoutes = () => (
  <>
    <Route index element={<Home />} />
    {routes.map(({ path, element: Element }) => (
      <Route key={path} path={path} element={<Element />} />
    ))}
    <Route path="*" element={<NotFound />} />
  </>
)
