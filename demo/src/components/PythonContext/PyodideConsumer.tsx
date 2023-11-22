import PyodideContext from './PyodideContext'
import { PyodideContextProps } from './PyodideContext'

export interface PyodideConsumerProps {
  children: (props: PyodideContextProps) => JSX.Element
}

const PyodideConsumer = ({ children }: PyodideConsumerProps) => (
  <PyodideContext.Consumer>
    {pyodideContextProps => children(pyodideContextProps)}
  </PyodideContext.Consumer>
)

export default PyodideConsumer
