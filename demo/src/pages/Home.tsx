import Box from '@mui/material/Box'
import Link from '@mui/material/Link'
import Typography from '@mui/material/Typography'

const Home = () => (
  <Box>
    <Box>
      <Typography variant="body2">
        This is a demonstration of option pricing formula implemented in Python.
      </Typography>
      <Typography variant="body2">
        You can find the project on github
        <Link href="https://github.com/rob-blackbourn/jetblack-options">
          {' here'}
        </Link>
        .
      </Typography>
    </Box>
  </Box>
)

export default Home
