import { ModelRoute } from './types'

import blackScholes73 from './blackScholes73'
import black76 from './black76'
import blackScholesMerton from './blackScholesMerton'
import garmanKohlhagen from './garmanKohlhagen'
import generalisedBlackScholes from './generalisedBlackScholes'
import baroneAdesiWhaley from './baroneAdesiWhaley'
import bjerksundStensland1993 from './bjerksundStensland1993'
import bjerksundStensland2002 from './bjerksundStensland2002'
import coxRossRubenstein from './coxRossRubenstein'
import europeanBinomialTree from './europeanBinomialTree'
import jarrowRudd from './jarrowRudd'
import leisenReimer from './leisenReimer'
import trinomialTree from './trinomialTree'

export const routes: ModelRoute[] = [
  blackScholes73,
  black76,
  blackScholesMerton,
  garmanKohlhagen,
  generalisedBlackScholes,
  baroneAdesiWhaley,
  bjerksundStensland1993,
  bjerksundStensland2002,
  coxRossRubenstein,
  europeanBinomialTree,
  jarrowRudd,
  leisenReimer,
  trinomialTree
]
