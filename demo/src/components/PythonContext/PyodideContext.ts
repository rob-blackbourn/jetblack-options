import React from 'react'

import { PyProxy, PyodideInterface } from 'pyodide'

export interface PyodideContextProps {
  pyodide: PyodideInterface | undefined
  micropip: PyProxy | undefined
  isRequirementsLoaded: boolean
}

const PyodideContext = React.createContext<PyodideContextProps>({
  pyodide: undefined,
  micropip: undefined,
  isRequirementsLoaded: false
})

export default PyodideContext
