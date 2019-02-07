#!/bin/bash

if [ "$#" -ne "2" ]; then
  echo "Usage: $0 <doSubmit> <do4params>"
  exit 0
fi

DATE=171201_unblind_datacards

DOSUBMIT=$1
DO4PARAMS=$2
infolder=output_run2ana_${DATE} #_datacard
outfolder=cards_run2ana_${DATE} #_datacard
do_tau_veto_unc=true
do_b_veto_unc=true
blind=false
#zvvstat=18
mkdir -p $outfolder

extraoptions="--do_ues=false" #--do_ggh=false --do_separate_qcdewk=false"

for channel in enu munu mumu ee nunu #taunu qcd
  do
  echo " ********************************"
  echo " *** Processing channel $channel"
  echo " ********************************"

  mkdir -p $outfolder/$channel

  echo "channel $channel"

  HistToIntegrate="jet2_pt"

  echo "channel $channel, HistToIntegrate: $HistToIntegrate"

  OUTNAME=$outfolder/$channel/vbfhinv_${channel}_13TeV.txt
  if (( "$DO4PARAMS" == "1" )); then
    OUTNAME=$outfolder/$channel/vbfhinv_${channel}_13TeV_4params.txt
  fi
  if (( "$DOSUBMIT" == "0" )); then
    echo "./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m 125 --channel $channel --do_latex true --do_datatop false --zvvstat 0 --qcdrate 0 --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --histoToIntegrate=$HistToIntegrate --do_b_veto_unc=$do_b_veto_unc --do_tau_veto_unc=$do_tau_veto_unc $extraoptions | tee $outfolder/$channel/card.log"
  else
    ./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m 125 --channel $channel --do_latex true --do_datatop false --zvvstat 0 --qcdrate 0 --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --histoToIntegrate=$HistToIntegrate --do_b_veto_unc=$do_b_veto_unc --do_tau_veto_unc=$do_tau_veto_unc $extraoptions | tee $outfolder/$channel/card_${channel}.log
  fi
done
