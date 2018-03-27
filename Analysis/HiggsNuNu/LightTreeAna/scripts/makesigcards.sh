#!/bin/bash

if [ "$#" -ne "2" ]; then
  echo "Usage: $0 <doSubmit> <do4params>"
  exit 0
fi

DATE=171201_unblind_datacards

DOSUBMIT=$1
DO4PARAMS=$2
infolder=output_run2ana_${DATE}
outfolder=cards_run2ana_${DATE}_masscards
do_tau_veto_unc=true
do_b_veto_unc=true
blind=false
zvvstat=0
mkdir -p $outfolder

extraoptions="--do_ues=false" #--do_ggh=false --do_separate_qcdewk=false"

for mass in 110 125 150 200 300 400 500 600 800 1000
  do
    mkdir -p $outfolder/$mass
    echo "------------------------"
    echo "---Processing mass $mass"
    echo "------------------------"
    OUTNAME=$outfolder/$mass/vbfhinv_${mass}_13TeV.txt
    if (( "$DO4PARAMS" == "1" )); then
	OUTNAME=$outfolder/$mass/vbfhinv_${mass}_13TeV_4params.txt
    fi
    if (( "$DOSUBMIT" == "0" )); then
	echo "./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m $mass --do_latex true --do_datatop false --zvvstat $zvvstat --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --histoToIntegrate=jet2_pt --do_b_veto_unc=$do_b_veto_unc --do_tau_veto_unc=$do_tau_veto_unc $extraoptions | tee $outfolder/$mass/card_${mass}.log"
    else
	./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m $mass --do_latex true --do_datatop false --zvvstat $zvvstat --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --histoToIntegrate=jet2_pt --do_b_veto_unc=$do_b_veto_unc --do_tau_veto_unc=$do_tau_veto_unc $extraoptions | tee $outfolder/$mass/card_${mass}.log
    fi
done

