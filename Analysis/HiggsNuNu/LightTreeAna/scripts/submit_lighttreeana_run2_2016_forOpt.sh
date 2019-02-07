#!/bin/sh
DOCERN=0
DOSUBMIT=1

DATE=171005_DetajjTrig_Opt

## Try and take the JOBWRAPPER and JOBSUBMIT commands
## from the environment if set, otherwise use these defaults
: ${JOBWRAPPER:="./scripts/generate_job.sh"}
: ${JOBSUBMIT:="eval"}

GRIDSETUP=1
if [ "$DOCERN" = "0" ]
  then
  JOBSCRIPT="./scripts/submit_ic_batch_job.sh"
else
  JOBSCRIPT="./scripts/submit_cern_batch_job.sh"
  GRIDSETUP=0
fi
export JOBSUBMIT=$JOBSCRIPT" "$JOBQUEUE


echo "Using job-wrapper: " $JOBWRAPPER
echo "Using job-submission: " $JOBSUBMIT

CONFIG=scripts/DefaultRun2Config_forOpt.cfg

QUEUEDIR=short #medium long

JOBDIRPREFIX=jobs_run2ana_${DATE}
JOBDIR=$JOBDIRPREFIX/
OUTPUTPREFIX=output_run2ana_${DATE}
OUTPUTDIR=$OUTPUTPREFIX/

OUTPUTNAME="output.root"

echo "Config file: $CONFIG"
mkdir -p $JOBDIR
mkdir -p $OUTPUTDIR

if [ "$DOCERN" = "0" ]
  then
  if [ "$QUEUEDIR" = "medium" ]
    then
    JOBQUEUE="5:59:0"
  elif [ "$QUEUEDIR" = "long" ]
    then
    JOBQUEUE="47:59:0"
  else
    JOBQUEUE="2:59:0"
  fi
else
  if [ "$QUEUEDIR" = "medium" ]
    then
    JOBQUEUE="1nd"
  elif [ "$QUEUEDIR" = "long" ]
    then
    JOBQUEUE="2nd"
  else
    JOBQUEUE="1nh"
  fi
fi
export JOBSUBMIT=$JOBSCRIPT" "$JOBQUEUE
echo "Using job-submission: " $JOBSUBMIT

echo "JOB name = $JOB"
for syst in "" JESUP JESDOWN TAUUP TAUDOWN BTAGUP BTAGDOWN LEPEFF_ELEUP LEPEFF_ELEDOWN LEPEFF_GSFUP LEPEFF_GSFDOWN LEPEFF_MUIDUP LEPEFF_MUIDDOWN LEPEFF_MUISOUP LEPEFF_MUISODOWN LEPEFF_MUTKUP LEPEFF_MUTKDOWN PUUP PUDOWN TRIGUP TRIGDOWN
#for syst in JESUP JESDOWN TAUUP TAUDOWN BTAGUP BTAGDOWN LEPEFF_ELEUP LEPEFF_ELEDOWN LEPEFF_GSFUP LEPEFF_GSFDOWN LEPEFF_MUIDUP LEPEFF_MUIDDOWN LEPEFF_MUISOUP LEPEFF_MUISODOWN LEPEFF_MUTKUP LEPEFF_MUTKDOWN PUUP PUDOWN TRIGUP TRIGDOWN
## Do not use till new JEC #NOTE TO RUN JER DOSMEAR MUST BE SET TO TRUE IN THE CONFIG
#for syst in JERBETTER JERWORSE
## For CMS_scale_met uncertainty
#for syst in UESUP UESDOWN
do
  mkdir -p $JOBDIR$syst
  mkdir -p $OUTPUTDIR$syst
  for channels in enu munu ee mumu nunu #taunu qcd
    do
    JOB=$channels
    #HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/weights.hists`
    #SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/weights.hists`
    #executable expect strings separated by "!"
    ## To produce all of the hist
    #HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/${channels}.hists`
    #SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/${channels}.hists`
    ## To produce all of the hist for datacard
    HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/${channels}_datacard.hists`
    SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/${channels}_datacard.hists`
    #HISTSTRING=`awk '{FS="\t"}{ORS="!"}{print $2}' scripts/${channels}_sig.hists`
    #SHAPESTRING=`awk '{ORS="!"}{print $1}' scripts/${channels}_sig.hists`
    ## To test for one hist
    #HISTSTRING=";#Delta#phi_{jj};Events"
    #SHAPESTRING="dijet_dphi(50,0.,3.1416)"
    #HISTSTRING=";E_{T,no-#mu}^{miss} (GeV);Events"
    #SHAPESTRING="metnomuons(25,200.,600.)"
    #HISTSTRING=";E_{T,no-#mu}^{miss} (GeV);Events!;Forward tag jet #eta;Events"
    #SHAPESTRING="metnomuons(25,200.,600.)!forward_tag_eta(25,-5.,5.)"
    #HISTSTRING=";#Delta#phi(E_{T,no-#mu}^{miss},j);Events"
    #SHAPESTRING="fourjetsmetnomu_mindphi(14,2.3,3.1416)"
    echo "Making histograms: " $SHAPESTRING
    OUTPUTNAME="$channels.root"
    MINDPHICUT="fourjetsmetnomu_mindphi\>=0.5"
    if [ "$channels" = "taunu" ]; then
	############MINDPHICUT="jetmetnomu_mindphi\>=1.0" #\&\&fourjetsmetnomu_mindphi\<2.3"
	#MINDPHICUT="jetmetnomu_mindphi\>=1.0\&\&fourjetsmetnomu_mindphi\<2.3"
      MINDPHICUT="fourjetsmetnomu_mindphi\>=0.5"
    fi
    if [ "$channels" = "qcd" ]; then
      MINDPHICUT="fourjetsmetnomu_mindphi\<0.5"
    fi
    if [ "$channels" = "ee" ]; then
      MINDPHICUT="fourjetsmetnoel_mindphi\>=0.5"
    fi
    if [ "$channels" = "enu" ]; then
      MINDPHICUT="fourjetsmetnoel_mindphi\>=0.5"
    fi
    if [ "$syst" = "" ]
      then
      $JOBWRAPPER "./bin/LTAnalysisRun2_2016 --cfg=$CONFIG --channel=$channels --histTitlePar='$HISTSTRING' --shapePar='$SHAPESTRING' -o $OUTPUTDIR$syst/$OUTPUTNAME --jetmetdphicut=$MINDPHICUT &> $JOBDIR$syst/$JOB.log" $JOBDIR$syst/$JOB.sh $GRIDSETUP
    else
	#if [ "$channels" = "nunu" ]
	   # then
	    #for other systs only run nunu
    #executable expect strings separated by "!"
	#HISTSTRING="Forward tag jet #eta;Events"
	#SHAPESTRING="forward_tag_eta(25,-5.,5.)"
	#echo "Making histograms: " $SHAPESTRING

      $JOBWRAPPER "./bin/LTAnalysisRun2_2016 --cfg=$CONFIG --channel=$channels --histTitlePar='$HISTSTRING' --shapePar='$SHAPESTRING' --syst=$syst -o $OUTPUTDIR$syst/$OUTPUTNAME --jetmetdphicut=$MINDPHICUT &> $JOBDIR$syst/$JOB.log" $JOBDIR$syst/$JOB.sh $GRIDSETUP
	#fi
    fi
    if [ "$DOSUBMIT" = "1" ]; then 
      $JOBSUBMIT $JOBDIR$syst/$JOB.sh
    else 
      echo "$JOBSUBMIT $JOBDIR$syst/$JOB.sh"
    fi
  done
done
