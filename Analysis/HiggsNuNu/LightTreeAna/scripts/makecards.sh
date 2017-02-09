#!/bin/bash

if [ "$#" -ne "2" ]; then
    echo "Usage: $0 <doSubmit> <do4params>"
    exit 0
fi


DOSUBMIT=$1
DO4PARAMS=$2
infolder=output_run2ana_170206_ICHEP/
outfolder=cards_run2ana_170206_ICHEP/
do_tau_veto_unc=false
blind=true
#zvvstat=18
mkdir -p $outfolder

extraoptions="--do_ues=false" #--do_ggh=false --do_separate_qcdewk=false"

for channel in qcd enu munu taunu mumu ee nunu #qcd #topl topb
do
    echo " ********************************"
    echo " *** Processing channel $channel"
    echo " ********************************"

    mkdir -p $outfolder/$channel

    echo "channel $channel"

    HistToIntegrate="jet2_pt"

    echo "channel $channel, HistToIntegrate: $HistToIntegrate"
	
    for minjjcut in 1101 #1601 #801 901 1001 1101 1201 1301 1401 1501 1601 1701 1801 1901
    do
	OUTNAME=$outfolder/$channel/vbfhinv_${channel}_13TeV.txt
	if (( "$DO4PARAMS" == "1" )); then
	    OUTNAME=$outfolder/$channel/vbfhinv_${channel}_13TeV_4params.txt
	fi
	if (( "$DOSUBMIT" == "0" )); then
	    echo "./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m 125 --channel $channel --do_latex true --do_datatop false --zvvstat 0 --qcdrate 0 --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --histoToIntegrate=$HistToIntegrate --do_tau_veto_unc=$do_tau_veto_unc $extraoptions | tee $outfolder/$channel/card.log"
	else
	    ./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m 125 --channel $channel --do_latex true --do_datatop false --zvvstat 0 --qcdrate 0 --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --histoToIntegrate=$HistToIntegrate --do_tau_veto_unc=$do_tau_veto_unc $extraoptions | tee $outfolder/$channel/card.log
	fi
    done
#done
done
