import CssBaseline from '@mui/material/CssBaseline'

import PythonApp from './components/PythonApp'

import { PyodideInterface } from 'pyodide'

declare global {
  interface Window {
    loadPyodide: () => Promise<PyodideInterface>
  }
}

function App() {
  return (
    <div>
      <CssBaseline />
      <PythonApp requirements={['jetblack-options']} />
    </div>
  )
}

export default App
