#!/bin/sh
command -v setupATLAS > /dev/null || (export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase && alias setupATLAS="source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh")
setupATLAS
lsetup git
asetup AthGeneration,21.6.36
export RIVET_ANALYSIS_PATH=$PWD
