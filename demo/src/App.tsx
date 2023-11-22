import Box from '@mui/material/Box'
import Container from '@mui/material/Container'
import CssBaseline from '@mui/material/CssBaseline'
import Link from '@mui/material/Link'
import Typography from '@mui/material/Typography'

import { PyodideProvider } from './components/PythonContext'

import { PyodideInterface } from 'pyodide'

import OptionValuer from './components/OptionValuer'

declare global {
  interface Window {
    loadPyodide: () => Promise<PyodideInterface>
  }
}

function App() {
  return (
    <div>
      <CssBaseline />
      <PyodideProvider requirements={['jetblack-options']}>
        <Container maxWidth="md" sx={{ width: '100%' }}>
          <Box>
            <Typography variant="h4">jetblack-options</Typography>
            <Typography variant="body1">
              This is a demonstration of option pricing formula implemented in
              Python.
            </Typography>
            <Link href="https://github.com/rob-blackbourn/jetblack-options">
              See here.
            </Link>
          </Box>
          <OptionValuer />
        </Container>
      </PyodideProvider>
    </div>
  )
}

export default App
