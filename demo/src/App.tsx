import React, { Suspense } from 'react'

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

import Layout from './pages/Layout'
import Home from './pages/Home'
import NotFound from './pages/NotFound'
import OptionRunner from './components/OptionRunner'

import { routes } from './nestedRoutes'

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
          <BrowserRouter basename="/jetblack-options">
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
          </BrowserRouter>
        </Container>
      </PyodideProvider>
    </ThemeProvider>
  )
}

export default App
