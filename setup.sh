#!/bin/sh
command -v setupATLAS > /dev/null || (export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase && alias setupATLAS="source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh")
setupATLAS
lsetup git
if [ -f .asetup.save ]; then
	asetup --restore
else
	asetup 21.6.67,AthGeneration
fi
export RIVET_ANALYSIS_PATH=$PWD
source setupRivet.sh
