# Useful scripts
This folder contains useful scripts to use on the FireSim platform.
The script buildAndInstall.sh should be on the following path ~/firesim/sw/firesim-software and it builds and install the workload which name is passed by parameter.
The script launchAndRun.sh should be on the following path ~/firesim/deploy and it executes the commands launchrunfarm, infrasetup and runworkload for the configuration file passed by parameter. If it doen't recieve any parameter, it executes the commands for the default config_runtime.ini configuration file.