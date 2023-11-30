import React from 'react'

import {
  Link as RouterLink,
  LinkProps as RouterLinkProps
} from 'react-router-dom'
import { LinkProps } from '@mui/material/Link'

import Container from '@mui/material/Container'
import CssBaseline from '@mui/material/CssBaseline'
import { ThemeProvider, createTheme } from '@mui/material/styles'

import { PyodideProvider } from './components/PythonContext'
import { PyodideInterface } from 'pyodide'

import RoutedApp from './RoutedApp'

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
          <RoutedApp />
        </Container>
      </PyodideProvider>
    </ThemeProvider>
  )
}

export default App
