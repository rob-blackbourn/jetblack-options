import React, { useEffect, useState } from 'react'

import { PyProxy, PyodideInterface } from 'pyodide'

import { Stack, CircularProgress } from '@mui/material'
import BlackScholes from './BlackScholes'

declare global {
  interface Window {
    loadPyodide: () => Promise<PyodideInterface>
  }
}

export interface PythonAppProps {
  requirements: string | string[]
}

const PythonApp: React.FC<PythonAppProps> = ({ requirements }) => {
  const [pyodide, setPyodide] = useState<PyodideInterface>()
  const [micropip, setMicropip] = useState<PyProxy>()
  const [isOptionsLoaded, setIsOptionsLoaded] = useState(false)

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

    micropip
      .install(requirements)
      .then(() => {
        console.log('Requirements installed')
        setIsOptionsLoaded(true)
      })
      .catch((error: Error) => {
        console.log(error)
      })
  }, [micropip, requirements])

  if (!(pyodide && micropip && isOptionsLoaded)) {
    return (
      <Stack>
        <CircularProgress />
      </Stack>
    )
  }

  console.log({ pyodide })

  return (
    <div>
      <BlackScholes pyodide={pyodide} />
    </div>
  )
}

export default PythonApp
