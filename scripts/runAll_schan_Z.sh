#!bin/sh

# colours
NC="\033[0m"
CYAN="\033[0;36m"
GREEN="\033[0;32m"

XSECFILE="xsec/xchan_Z.txt"
# remove cross section file from previous runs
[ -f $XSECFILE ] && rm -f $XSECFILE

MASSES=(100 200 300 400 500 600 700 800)

BASECMD="rivet --analysis=SSHiggs --nominal-weight Weight_MERGING=0.000 --skip-weights"
OUTPUTPREFIX="outputs/schan_ZGamma_ll_eemumu_MHPPL_"
INPUTFILEBASE="/eos/user/s/sdysch/SSHiggs/outputs/schan_ZGamma_ll_eemumu_MHPPL_"

# option to debug with $1 events
[ $1 ] && INFO "Restricting to only $1 events" && BASECMD="$BASECMD -n $1"

for mass in ${MASSES[@]}; do
	OUTPUTFILE="${OUTPUTPREFIX}${mass}GEV.yoda"
	INPUTFILE="${INPUTFILEBASE}${mass}GEV/Events/run_01/tag_1_pythia8_events.hepmc.gz"
	SUFFIX="--histo-file $OUTPUTFILE --pwd ${INPUTFILE}"

	CMD="${BASECMD} ${SUFFIX}"
	#echo $CMD

	LOGFILE="/tmp/sdysch/log_${mass}"
	[ -f $LOGFILE ] && rm -f $LOGFILE
	$CMD | tee $LOGFILE
	#$CMD

	# parse log file for cross section
	grep "Rivet cross section:" $LOGFILE | awk -F ": " '{print $2}' >> $XSECFILE
done

function INFO() {
	echo -e "${GREEN}[INFO] $1 ${NC}"
}
