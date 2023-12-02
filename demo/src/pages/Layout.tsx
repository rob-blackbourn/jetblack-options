import { useState, useEffect } from 'react'
import { Outlet, useNavigate } from 'react-router-dom'

import Box from '@mui/material/Box'
import FormControl from '@mui/material/FormControl'
import InputLabel from '@mui/material/InputLabel'
import MenuItem from '@mui/material/MenuItem'
import Select, { SelectChangeEvent } from '@mui/material/Select'
import Stack from '@mui/material/Stack'
import Typography from '@mui/material/Typography'

import { routes } from '../nestedRoutes'

const Layout = () => {
  const [page, setPage] = useState<string>('')
  const navigate = useNavigate()

  const handlePageChange = (event: SelectChangeEvent) => {
    setPage(event.target.value as string)
  }

  useEffect(() => {
    if (!page) {
      return
    }

    navigate(page)
  }, [page, navigate])

  return (
    <Stack>
      <Box sx={{ minWidth: 120, mt: 2, mb: 2 }}>
        <Typography variant="h4">jetblack-options</Typography>

        <FormControl fullWidth sx={{ mt: 2 }}>
          <InputLabel id="page-select-label">Page</InputLabel>
          <Select
            labelId="page-select-label"
            id="page-select"
            value={page}
            label="Page"
            onChange={handlePageChange}
          >
            {routes.map(({ path, schema: { title } }) => (
              <MenuItem key={path} value={path}>
                {title}
              </MenuItem>
            ))}
          </Select>
        </FormControl>
      </Box>

      <Outlet />
    </Stack>
  )
}

export default Layout
