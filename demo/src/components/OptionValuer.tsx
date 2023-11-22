import React, { useState } from 'react'

import Box from '@mui/material/Box'
import Container from '@mui/material/Container'
import Tabs from '@mui/material/Tabs'
import Tab from '@mui/material/Tab'

import BlackScholes73 from './BlackScholes73'
import Black76 from './Black76'
import BlackScholesMerton from './BlackScholesMerton'
import GeneralisedBlackScholes from './GeneralisedBlackScholes'
import GarmanKolhagen from './GarmanKolhagen'

interface TabPanelProps {
  children?: React.ReactNode
  index: number
  value: number
}

const TabPanel: React.FC<TabPanelProps> = ({
  children,
  value,
  index,
  ...other
}) => (
  <div
    role="tabpanel"
    hidden={value !== index}
    id={`simple-tabpanel-${index}`}
    aria-labelledby={`simple-tab-${index}`}
    {...other}
  >
    {value === index && <Box sx={{ p: 3 }}>{children}</Box>}
  </div>
)

function a11yProps(index: number) {
  return {
    id: `simple-tab-${index}`,
    'aria-controls': `simple-tabpanel-${index}`
  }
}

export interface OptionValuerProps {}

const PythonApp: React.FC<OptionValuerProps> = () => {
  const [tabIndex, setTabIndex] = useState(0)

  const handleTabIndexChange = (
    _event: React.SyntheticEvent,
    newValue: number
  ) => {
    setTabIndex(newValue)
  }

  return (
    <Container maxWidth="md" sx={{ width: '100%' }}>
      <Box sx={{ borderBottom: 1, borderColor: 'divider' }}>
        <Tabs value={tabIndex} onChange={handleTabIndexChange}>
          <Tab label="Black-Scholes 73" {...a11yProps(0)} />
          <Tab label="Black 76" {...a11yProps(1)} />
          <Tab label="Black-Scholes-Merton" {...a11yProps(2)} />
          <Tab label="Generalised Black-Scholes" {...a11yProps(3)} />
          <Tab label="Garman Kolhagen" {...a11yProps(4)} />
        </Tabs>
      </Box>
      <TabPanel value={tabIndex} index={0}>
        <BlackScholes73 />
      </TabPanel>
      <TabPanel value={tabIndex} index={1}>
        <Black76 />
      </TabPanel>
      <TabPanel value={tabIndex} index={2}>
        <BlackScholesMerton />
      </TabPanel>
      <TabPanel value={tabIndex} index={3}>
        <GeneralisedBlackScholes />
      </TabPanel>
      <TabPanel value={tabIndex} index={4}>
        <GarmanKolhagen />
      </TabPanel>
    </Container>
  )
}

export default PythonApp
