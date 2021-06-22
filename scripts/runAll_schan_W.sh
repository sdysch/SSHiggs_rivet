#!bin/sh

# colours
NC="\033[0m"
CYAN="\033[0;36m"
GREEN="\033[0;32m"


BASECMD="rivet --analysis=SSHiggs --nominal-weight Weight_MERGING=0.000 --skip-weights"
INPUTFILEBASE="/eos/user/s/sdysch/SSHiggs/outputs/schan_W_ll_eemumu/Events/"

# option to debug with $1 events
[ $1 ] && INFO "Restricting to only $1 events" && BASECMD="$BASECMD -n $1"

CMD100GEV="${BASECMD} --histo-file outputs/schan_W_ll_eemumu_MHPPL_100GEV.yoda --pwd ${INPUTFILEBASE}/run_01/tag_1_pythia8_events.hepmc.gz"
CMD200GEV="${BASECMD} --histo-file outputs/schan_W_ll_eemumu_MHPPL_200GEV.yoda --pwd ${INPUTFILEBASE}/run_02/tag_1_pythia8_events.hepmc.gz"
CMD300GEV="${BASECMD} --histo-file outputs/schan_W_ll_eemumu_MHPPL_300GEV.yoda --pwd ${INPUTFILEBASE}/run_03/tag_1_pythia8_events.hepmc.gz"
CMD400GEV="${BASECMD} --histo-file outputs/schan_W_ll_eemumu_MHPPL_400GEV.yoda --pwd ${INPUTFILEBASE}/run_04/tag_1_pythia8_events.hepmc.gz"
CMD500GEV="${BASECMD} --histo-file outputs/schan_W_ll_eemumu_MHPPL_500GEV.yoda --pwd ${INPUTFILEBASE}/run_05/tag_1_pythia8_events.hepmc.gz"
CMD600GEV="${BASECMD} --histo-file outputs/schan_W_ll_eemumu_MHPPL_600GEV.yoda --pwd ${INPUTFILEBASE}/run_06/tag_1_pythia8_events.hepmc.gz"
CMD700GEV="${BASECMD} --histo-file outputs/schan_W_ll_eemumu_MHPPL_700GEV.yoda --pwd ${INPUTFILEBASE}/run_07/tag_1_pythia8_events.hepmc.gz"
CMD800GEV="${BASECMD} --histo-file outputs/schan_W_ll_eemumu_MHPPL_800GEV.yoda --pwd ${INPUTFILEBASE}/run_08/tag_1_pythia8_events.hepmc.gz"

$CMD100GEV
$CMD200GEV
$CMD300GEV
$CMD400GEV
$CMD500GEV
$CMD600GEV
$CMD700GEV
$CMD800GEV

function INFO() {
	echo -e "${GREEN}[INFO] $1 ${NC}"
}
