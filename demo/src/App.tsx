import React, { Suspense, lazy } from 'react'

import {
  BrowserRouter,
  Route,
  Link as RouterLink,
  LinkProps as RouterLinkProps,
  Routes
} from 'react-router-dom'
import { LinkProps } from '@mui/material/Link'

import Container from '@mui/material/Container'
import CssBaseline from '@mui/material/CssBaseline'
import { ThemeProvider, createTheme } from '@mui/material/styles'

import { PyodideProvider } from './components/PythonContext'
import Loading from './components/Loading'
import { PyodideInterface } from 'pyodide'

import Home from './pages/Home'
import Layout from './pages/Layout'
import NotFound from './pages/NotFound'

const Black76 = lazy(() => import('./components/Black76'))
const BlackScholes73 = lazy(() => import('./components/BlackScholes73'))
const BlackScholesMerton = lazy(() => import('./components/BlackScholesMerton'))
const GarmanKohlhagen = lazy(() => import('./components/GarmanKohlhagen'))

declare global {
  interface Window {
    loadPyodide: () => Promise<PyodideInterface>
  }
}

const LinkBehavior = React.forwardRef<
  HTMLAnchorElement,
  Omit<RouterLinkProps, 'to'> & { href: RouterLinkProps['to'] }
>((props, ref) => {
  const { href, ...other } = props
  // Map href (Material UI) -> to (react-router)
  return <RouterLink ref={ref} to={href} {...other} />
})

const theme = createTheme({
  components: {
    MuiLink: {
      defaultProps: {
        component: LinkBehavior
      } as LinkProps
    },
    MuiButtonBase: {
      defaultProps: {
        LinkComponent: LinkBehavior
      }
    }
  }
})

function App() {
  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <PyodideProvider requirements={['jetblack-options']}>
        <Container maxWidth="md" sx={{ width: '100%' }}>
          <BrowserRouter>
            <Suspense fallback={<Loading />}>
              <Routes>
                <Route path="/" element={<Layout />}>
                  <Route index element={<Home />} />
                  <Route path="/black-76" element={<Black76 />} />
                  <Route
                    path="/black-scholes-73"
                    element={<BlackScholes73 />}
                  />
                  <Route
                    path="/black-scholes-merton"
                    element={<BlackScholesMerton />}
                  />
                  <Route
                    path="/garman-kohlhagen"
                    element={<GarmanKohlhagen />}
                  />
                  <Route path="*" element={<NotFound />} />
                </Route>
              </Routes>
            </Suspense>
          </BrowserRouter>
        </Container>
      </PyodideProvider>
    </ThemeProvider>
  )
}

export default App
