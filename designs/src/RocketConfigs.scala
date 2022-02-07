package chipyard

import freechips.rocketchip.config.{Config}
import freechips.rocketchip.diplomacy.{AsynchronousCrossing}

// --------------
// Rocket Configs
// --------------

class RocketConfig extends Config(
  new freechips.rocketchip.subsystem.WithNBigCores(1) ++         
  new chipyard.config.AbstractConfig)

class HwachaRocketConfig extends Config(
  new chipyard.config.WithHwachaTest ++
  new hwacha.DefaultHwachaConfig ++                              
  new freechips.rocketchip.subsystem.WithNBigCores(1) ++
  new chipyard.config.AbstractConfig)

class GemminiRocketConfig extends Config(
  new gemmini.DefaultGemminiConfig ++                            
  new freechips.rocketchip.subsystem.WithNBigCores(1) ++
  new chipyard.config.AbstractConfig)

class QuadRocketConfig extends Config(
  new freechips.rocketchip.subsystem.WithNBigCores(4) ++    
  new chipyard.config.AbstractConfig)

class ThreeRocketConfig extends Config(
  new freechips.rocketchip.subsystem.WithNBigCores(3) ++    
  new chipyard.config.AbstractConfig)

class DualRocketConfig extends Config(
  new freechips.rocketchip.subsystem.WithNBigCores(2) ++    
  new chipyard.config.AbstractConfig)

