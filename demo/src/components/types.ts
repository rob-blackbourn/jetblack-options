export interface Greeks {
  delta?: number
  gamma?: number
  theta?: number
  vega?: number
  rho?: number
}

export interface OptionResults {
  price: number
  analytic: Greeks
  numeric: Greeks
}
