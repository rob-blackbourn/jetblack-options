import { PyodideInterface } from 'pyodide'
import { OptionResults } from './types'

export function runBlack76(
  pyodide: PyodideInterface,
  args: object,
  analyticImportPath: string,
  numericImportPath: string,
  priceArgNames: string[],
  analyticGreeks: Record<string, string[] | null>
): OptionResults {
  const locals = pyodide.toPy({ args })

  const analyticImports = ['price', ...Object.keys(analyticGreeks)].join(', ')

  const analyticValuations = `{ ${Object.entries(analyticGreeks)
    .map(([greek, args]) =>
      args == null
        ? `${greek}: None`
        : `'${greek}': ${greek}(${args?.join(', ')})`
    )
    .join(', ')} }`

  const extractArgs = Object.keys(args)
    .map(name => `${name} = args['${name}']`)
    .join('\n')

  const script = `
from ${analyticImportPath} import (${analyticImports})
from ${numericImportPath} import NumericGreeks

${extractArgs}

analytics = ${analyticValuations}

ng = NumericGreeks(price)

numerics = { ${Object.keys(analyticGreeks).map(
    greek => `'${greek}': ng.${greek}(${priceArgNames.join(', ')})`
  )} }

{
    'price': price(${priceArgNames.join(', ')}),
    'analytic': analytics,
    'numeric': numerics,
}`
  console.log(script)

  const results = pyodide.runPython(script, { locals })

  const dct = results.toJs()
  const analytic = dct.get('analytic')
  const numeric = dct.get('numeric')
  const optionResults: OptionResults = {
    price: dct.get('price'),
    analytic: Object.keys(analyticGreeks).reduce(
      (obj, key) => ({ ...obj, [key]: analytic.get(key) }),
      {}
    ),
    numeric: Object.keys(analyticGreeks).reduce(
      (obj, key) => ({ ...obj, [key]: numeric.get(key) }),
      {}
    )
  }

  results.destroy()

  return optionResults
}
