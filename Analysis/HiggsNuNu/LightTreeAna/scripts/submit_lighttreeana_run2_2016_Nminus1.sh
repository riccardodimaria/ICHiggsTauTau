#!/bin/sh
DOSUBMIT=1

#DATE=171215_unblind_Nminus1_plots
DATE=180207_unblind_ARC_Nminus1_plots

## Try and take the JOBWRAPPER and JOBSUBMIT commands
## from the environment if set, otherwise use these defaults
: ${JOBWRAPPER:="./scripts/generate_job.sh"}
: ${JOBSUBMIT:="eval"}

GRIDSETUP=1
JOBSCRIPT="./scripts/submit_ic_batch_job.sh"

export JOBSUBMIT=$JOBSCRIPT" "$JOBQUEUE


echo "Using job-wrapper: " $JOBWRAPPER
echo "Using job-submission: " $JOBSUBMIT

#CONFIG=scripts/DefaultRun2Config_Nminus1_met.cfg
CONFIG=scripts/DefaultRun2Config_Nminus1_dphi.cfg
#CONFIG=scripts/DefaultRun2Config_Nminus1_deta.cfg
#CONFIG=scripts/DefaultRun2Config_Nminus1_mindphi.cfg
#CONFIG=scripts/DefaultRun2Config_Nminus1_dijet.cfg

QUEUEDIR=short #medium long
# no_METno_cut no_dijet_dphi_cut no_dijet_deta_cut no_fourjetsmetnomu_mindphi_cut no_dijet_M_cut
JOBDIRPREFIX=jobs_run2ana_${DATE}/
JOBDIR=$JOBDIRPREFIX/
OUTPUTPREFIX=output_run2ana_${DATE}/
OUTPUTDIR=$OUTPUTPREFIX/

OUTPUTNAME="output.root"

echo "Config file: $CONFIG"
mkdir -p $JOBDIR
mkdir -p $OUTPUTDIR

if [ "$QUEUEDIR" = "medium" ]
  then
  JOBQUEUE="5:59:0"
elif [ "$QUEUEDIR" = "long" ]
  then
  JOBQUEUE="47:59:0"
else
  JOBQUEUE="2:59:0"
fi

export JOBSUBMIT=$JOBSCRIPT" "$JOBQUEUE
echo "Using job-submission: " $JOBSUBMIT

echo "JOB name = $JOB"
for syst in ""
do
  mkdir -p $JOBDIR$syst
  mkdir -p $OUTPUTDIR$syst
  for channels in nunu ee mumu enu munu
    do
    JOB=$channels
    #executable expect strings separated by "!"

    ## no_METno_cut
#     HISTSTRING=";E_{T,no-#mu}^{miss} [GeV];Events"
#     SHAPESTRING="metnomuons(50,0.,1000.)"

    ## no_dijet_dphi_cut
    HISTSTRING=";#Delta#phi_{jj};Events"
    SHAPESTRING="dijet_dphi(30,0.,3.1416)"
    #SHAPESTRING="dijet_dphi(60,0.,3.1416)"

    ## no_dijet_deta_cut
#     HISTSTRING=";#Delta#eta_{jj};Events"
#     SHAPESTRING="dijet_deta(35,1.,8.)"

    ## no_fourjetsmetnomu_mindphi_cut
#     HISTSTRING=";#Delta#phi(E_{T,no-#mu}^{miss},j);Events"
#     SHAPESTRING="fourjetsmetnomu_mindphi(24,0.,3.1416)"

    ## no_dijet_M_cut
#     HISTSTRING=";M_{jj} [GeV];Events"
#     SHAPESTRING="dijet_M(50,300.,5300.)"

    echo "Making histograms: " $SHAPESTRING
    OUTPUTNAME="$channels.root"
    MINDPHICUT="fourjetsmetnomu_mindphi\>=0.5"

    if [ "$syst" = "" ]
      then
      $JOBWRAPPER "./bin/LTAnalysisRun2_2016 --cfg=$CONFIG --channel=$channels --histTitlePar='$HISTSTRING' --shapePar='$SHAPESTRING' -o $OUTPUTDIR$syst/$OUTPUTNAME --jetmetdphicut=$MINDPHICUT &> $JOBDIR$syst/$JOB.log" $JOBDIR$syst/$JOB.sh $GRIDSETUP
    else

      $JOBWRAPPER "./bin/LTAnalysisRun2_2016 --cfg=$CONFIG --channel=$channels --histTitlePar='$HISTSTRING' --shapePar='$SHAPESTRING' --syst=$syst -o $OUTPUTDIR$syst/$OUTPUTNAME --jetmetdphicut=$MINDPHICUT &> $JOBDIR$syst/$JOB.log" $JOBDIR$syst/$JOB.sh $GRIDSETUP

    fi
    if [ "$DOSUBMIT" = "1" ]; then 
      $JOBSUBMIT $JOBDIR$syst/$JOB.sh
    else 
      echo "$JOBSUBMIT $JOBDIR$syst/$JOB.sh"
    fi
  done
done
