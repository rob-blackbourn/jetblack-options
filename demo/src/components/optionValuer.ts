import { PyodideInterface } from 'pyodide'
import { OptionResults } from './types'

export function valueOption(
  pyodide: PyodideInterface,
  args: object,
  analyticImportPath: string,
  pricePrototype: string[],
  analyticsPrototypes: Record<string, string[] | null>,
  bumpFactoryPrototype: string[],
  bumpPrototype: string[]
): OptionResults {
  const locals = pyodide.toPy({ args })

  const analyticImports = [
    'price',
    'make_bumper',
    ...Object.entries(analyticsPrototypes)
      .filter(([, args]) => args != null)
      .map(([name]) => name)
  ].join(', ')

  const analyticValuations = `{ ${Object.entries(analyticsPrototypes)
    .map(([greek, args]) =>
      args == null
        ? `'${greek}': None`
        : `'${greek}': ${greek}(${args?.join(', ')})`
    )
    .join(', ')} }`

  const extractArgs = Object.keys(args)
    .map(name => `${name} = args['${name}']`)
    .join('\n')

  const script = `
from ${analyticImportPath} import (${analyticImports})

${extractArgs}

analytics = ${analyticValuations}

ng = make_bumper(${bumpFactoryPrototype.join(', ')})

numerics = { ${Object.keys(analyticsPrototypes).map(
    greek => `'${greek}': ng.${greek}(${bumpPrototype.join(', ')})`
  )} }

{
    'price': price(${pricePrototype.join(', ')}),
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
    analytic: Object.keys(analyticsPrototypes).reduce(
      (obj, key) => ({ ...obj, [key]: analytic.get(key) }),
      {}
    ),
    numeric: Object.keys(analyticsPrototypes).reduce(
      (obj, key) => ({ ...obj, [key]: numeric.get(key) }),
      {}
    )
  }

  results.destroy()

  return optionResults
}
