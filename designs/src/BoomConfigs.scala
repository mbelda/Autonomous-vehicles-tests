package chipyard

import freechips.rocketchip.config.{Config}

// ---------------------
// BOOM Configs
// ---------------------

class LargeBoomConfig extends Config(
  new boom.common.WithNLargeBooms(1) ++                          
  new chipyard.config.AbstractConfig)

class HwachaLargeBoomConfig extends Config(
  new chipyard.config.WithHwachaTest ++
  new hwacha.DefaultHwachaConfig ++                              
  new boom.common.WithNLargeBooms(1) ++
  new chipyard.config.AbstractConfig)

class ThreeLargeBoomConfig extends Config(
  new boom.common.WithNLargeBooms(3) ++                          
  new chipyard.config.AbstractConfig)

class DualBoomConfig extends Config(
  new boom.common.WithNLargeBooms(2) ++                          
  new chipyard.config.AbstractConfig)

class GemminiBoomConfig extends Config(
  new gemmini.DefaultGemminiConfig ++                            
  new boom.common.WithNLargeBooms(1) ++
  new chipyard.config.AbstractConfig)
