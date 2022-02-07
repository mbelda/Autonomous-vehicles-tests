#!/bin/bash
if [ $# -eq 0 ]
  then
    firesim launchrunfarm
    firesim infrasetup
    firesim runworkload
else
    firesim launchrunfarm -c $1
    firesim infrasetup -c $1
    firesim runworkload -c $1
fi
