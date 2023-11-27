import { ModelRoute } from './types'

import blackScholes73 from './blackScholes73'
import black76 from './black76'
import blackScholesMerton from './blackScholesMerton'
import garmanKohlhagen from './garmanKohlhagen'
import generalisedBlackScholes from './generalisedBlackScholes'
import baroneAdesiWhaley from './baroneAdesiWhaley'

export const routes: ModelRoute[] = [
  blackScholes73,
  black76,
  blackScholesMerton,
  garmanKohlhagen,
  generalisedBlackScholes,
  baroneAdesiWhaley
]
