package chipyard

import freechips.rocketchip.config.{Config}

// ---------------------
// Heterogenous Configs
// ---------------------


class HwachaLargeBoomAndHwachaRocketConfig extends Config(
  new chipyard.config.WithHwachaTest ++
  new hwacha.DefaultHwachaConfig ++                          
  new boom.common.WithNLargeBooms(1) ++                      
  new freechips.rocketchip.subsystem.WithNBigCores(1) ++     
  new chipyard.config.AbstractConfig)
