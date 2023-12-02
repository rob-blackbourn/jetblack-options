import React from 'react'

import Paper from '@mui/material/Paper'
import Table from '@mui/material/Table'
import TableBody from '@mui/material/TableBody'
import TableCell from '@mui/material/TableCell'
import TableContainer from '@mui/material/TableContainer'
import TableHead from '@mui/material/TableHead'
import TableRow from '@mui/material/TableRow'
import { OptionResults } from './types'
import { SxProps, Theme } from '@mui/material'

const formatNumber = new Intl.NumberFormat(undefined, {
  minimumFractionDigits: 5,
  maximumFractionDigits: 5
}).format

const formatPercent = new Intl.NumberFormat(undefined, {
  style: 'percent',
  minimumFractionDigits: 5,
  maximumFractionDigits: 5
}).format

const formatOptionalNumber = (value: number | undefined) =>
  value == null ? '' : formatNumber(value)
const formatOptionalPercent = (value: number | undefined) =>
  value == null ? '' : formatPercent(value)
const calcDiff = (
  lhs: number | undefined,
  rhs: number | undefined
): number | undefined =>
  lhs === undefined || rhs === undefined ? undefined : rhs - lhs
const calcPercentageDiff = (
  lhs: number | undefined,
  rhs: number | undefined
): number | undefined =>
  lhs === undefined || rhs === undefined ? undefined : (rhs - lhs) / lhs

export interface OptionResultViewProps {
  optionResults: OptionResults
  sx?: SxProps<Theme>
}

const OptionResultView: React.FC<OptionResultViewProps> = ({
  optionResults,
  sx = {}
}) => {
  return (
    <TableContainer component={Paper} sx={sx}>
      <Table sx={{ minWidth: 650 }} aria-label="simple table">
        <TableHead>
          <TableRow>
            <TableCell>Name</TableCell>
            <TableCell align="right">Analytic</TableCell>
            <TableCell align="right">Numeric</TableCell>
            <TableCell align="right">diff</TableCell>
            <TableCell align="right">% diff</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          <TableRow>
            <TableCell component="th" scope="row">
              Price
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.price)}
            </TableCell>
            <TableCell align="right"></TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Delta
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.analytic.delta)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.numeric.delta)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(
                calcDiff(
                  optionResults.analytic.delta,
                  optionResults.numeric.delta
                )
              )}
            </TableCell>
            <TableCell align="right">
              {formatOptionalPercent(
                calcPercentageDiff(
                  optionResults.numeric.delta,
                  optionResults.analytic.delta
                )
              )}
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Gamma
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.analytic.gamma)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.numeric.gamma)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(
                calcDiff(
                  optionResults.analytic.gamma,
                  optionResults.numeric.gamma
                )
              )}
            </TableCell>
            <TableCell align="right">
              {formatOptionalPercent(
                calcPercentageDiff(
                  optionResults.numeric.gamma,
                  optionResults.analytic.gamma
                )
              )}
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Theta
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.analytic.theta)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.numeric.theta)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(
                calcDiff(
                  optionResults.analytic.theta,
                  optionResults.numeric.theta
                )
              )}
            </TableCell>
            <TableCell align="right">
              {formatOptionalPercent(
                calcPercentageDiff(
                  optionResults.numeric.theta,
                  optionResults.analytic.theta
                )
              )}
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Vega
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.analytic.vega)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.numeric.vega)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(
                calcDiff(
                  optionResults.analytic.vega,
                  optionResults.numeric.vega
                )
              )}
            </TableCell>
            <TableCell align="right">
              {formatOptionalPercent(
                calcPercentageDiff(
                  optionResults.numeric.vega,
                  optionResults.analytic.vega
                )
              )}
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Rho
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.analytic.rho)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(optionResults.numeric.rho)}
            </TableCell>
            <TableCell align="right">
              {formatOptionalNumber(
                calcDiff(optionResults.analytic.rho, optionResults.numeric.rho)
              )}
            </TableCell>
            <TableCell align="right">
              {formatOptionalPercent(
                calcPercentageDiff(
                  optionResults.numeric.rho,
                  optionResults.analytic.rho
                )
              )}
            </TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>
  )
}

export default OptionResultView
