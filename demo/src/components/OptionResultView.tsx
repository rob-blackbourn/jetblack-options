import React from 'react'

import Paper from '@mui/material/Paper'
import Table from '@mui/material/Table'
import TableBody from '@mui/material/TableBody'
import TableCell from '@mui/material/TableCell'
import TableContainer from '@mui/material/TableContainer'
import TableHead from '@mui/material/TableHead'
import TableRow from '@mui/material/TableRow'
import { OptionResults } from './types'

const formatNumber = new Intl.NumberFormat(undefined, {
  minimumFractionDigits: 5,
  maximumFractionDigits: 5
}).format

const formatPercent = new Intl.NumberFormat(undefined, {
  style: 'percent',
  minimumFractionDigits: 5,
  maximumFractionDigits: 5
}).format

export interface OptionResultViewProps {
  optionResults: OptionResults
}

const OptionResultView: React.FC<OptionResultViewProps> = ({
  optionResults
}) => {
  return (
    <TableContainer component={Paper}>
      <Table sx={{ minWidth: 650 }} aria-label="simple table" size="small">
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
              {formatNumber(optionResults.price)}
            </TableCell>
            <TableCell align="right"></TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Delta
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.analytic.delta)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.numeric.delta)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(
                optionResults.numeric.delta - optionResults.numeric.delta
              )}
            </TableCell>
            <TableCell align="right">
              {formatPercent(
                (optionResults.numeric.delta - optionResults.analytic.delta) /
                  optionResults.analytic.delta
              )}
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Gamma
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.analytic.gamma)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.numeric.gamma)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(
                optionResults.analytic.gamma - optionResults.numeric.gamma
              )}
            </TableCell>
            <TableCell align="right">
              {formatPercent(
                (optionResults.numeric.gamma - optionResults.analytic.gamma) /
                  optionResults.analytic.gamma
              )}
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Theta
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.analytic.theta)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.numeric.theta)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(
                optionResults.analytic.theta - optionResults.numeric.theta
              )}
            </TableCell>
            <TableCell align="right">
              {formatPercent(
                (optionResults.numeric.theta - optionResults.analytic.theta) /
                  optionResults.analytic.theta
              )}
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Vega
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.analytic.vega)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.numeric.vega)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(
                optionResults.analytic.vega - optionResults.numeric.vega
              )}
            </TableCell>
            <TableCell align="right">
              {formatPercent(
                (optionResults.numeric.vega - optionResults.analytic.vega) /
                  optionResults.analytic.vega
              )}
            </TableCell>
          </TableRow>
          <TableRow>
            <TableCell component="th" scope="row">
              Rho
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.analytic.rho)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(optionResults.numeric.rho)}
            </TableCell>
            <TableCell align="right">
              {formatNumber(
                optionResults.analytic.rho - optionResults.numeric.rho
              )}
            </TableCell>
            <TableCell align="right">
              {formatPercent(
                (optionResults.numeric.rho - optionResults.analytic.rho) /
                  optionResults.analytic.rho
              )}
            </TableCell>
          </TableRow>
        </TableBody>
      </Table>
    </TableContainer>
  )
}

export default OptionResultView
