import React, { useEffect, useState } from 'react'

import PyodideContext from './PyodideContext'
import { PyProxy, PyodideInterface } from 'pyodide'

declare global {
  interface Window {
    loadPyodide: () => Promise<PyodideInterface>
  }
}

export interface PyodideProviderProps {
  requirements?: string | string[]
  children?: React.ReactNode
}

const PyodideProvider: React.FC<PyodideProviderProps> = ({
  requirements,
  children
}) => {
  const [pyodide, setPyodide] = useState<PyodideInterface>()
  const [micropip, setMicropip] = useState<PyProxy>()
  const [isRequirementsLoaded, setIsRequirementsLoaded] = useState(false)

  useEffect(() => {
    window
      .loadPyodide()
      .then(pyodide => setPyodide(pyodide))
      .catch(error => console.log(error))
  }, [])

  useEffect(() => {
    if (!pyodide) {
      return
    }

    pyodide
      .loadPackage('micropip')
      .then(() => {
        console.log('Micropip loaded')
        const micropip = pyodide.pyimport('micropip')
        setMicropip(micropip)
      })
      .catch(error => {
        console.log(error)
      })
  }, [pyodide])

  useEffect(() => {
    if (!micropip) {
      return
    }

    if (!requirements) {
      setIsRequirementsLoaded(true)
      return
    }

    micropip
      .install(requirements)
      .then(() => {
        console.log('Requirements installed')
        setIsRequirementsLoaded(true)
      })
      .catch((error: Error) => {
        console.log(error)
      })
  }, [micropip, requirements])

  return (
    <PyodideContext.Provider
      value={{
        pyodide,
        micropip,
        isRequirementsLoaded
      }}
    >
      {children}
    </PyodideContext.Provider>
  )
}

export default PyodideProvider
