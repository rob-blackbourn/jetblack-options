import CssBaseline from '@mui/material/CssBaseline'

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
        <OptionValuer />
      </PyodideProvider>
    </div>
  )
}

export default App
