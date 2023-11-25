import Stack from '@mui/material/Stack'
import Typography from '@mui/material/Typography'

import WarningAmberIcon from '@mui/icons-material/WarningAmber'

const NotFound = () => (
  <Stack direction="row">
    <WarningAmberIcon color="warning" fontSize="large" />
    <Typography variant="h4">File not found</Typography>
  </Stack>
)

export default NotFound
