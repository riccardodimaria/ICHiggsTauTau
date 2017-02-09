#!/bin/sh
DOCERN=0
DOSUBMIT=0
#JETTYPE="ak4SlimmedJetsPuppi"
JETTYPE="pfJetsPFlow"
MYEXEC=LightTreeMakerFromMiniAODRDM
PRODUCTION=161031
PRODUSER=rdimaria
JPTCUTVAL=40
## Try and take the JOBWRAPPER and JOBSUBMIT commands
## from the environment if set, otherwise use these defaults
: ${JOBWRAPPER:="./scripts/generate_job.sh $DOCERN $MYEXEC $PRODUCTION"}
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

INPUTPARAMS="filelists/$PRODUCTION/Params${PRODUCTION}.dat"
CONFIG=scripts/DefaultLightTreeConfig_mc.cfg


for SYST in central #JESUP JESDOWN JERBETTER JERWORSE UESUP UESDOWN ELEEFFUP ELEEFFDOWN MUEFFUP MUEFFDOWN #NOTE TO RUN JER DOSMEAR MUST BE SET TO TRUE IN THE CONFIG
#for SYST in JESUP JESDOWN JERBETTER
#for SYST in JERWORSE UESUP UESDOWN
#for SYST in ELEEFFUP ELEEFFDOWN
#for SYST in MUEFFUP MUEFFDOWN
  do
  SYSTOPTIONS="--dojessyst=false --dojersyst=false"

  #JOBDIRPREFIX=/vols/cms/rd1715/HiggsToInv/jobs_lighttree_${PRODUCTION}_ICHEP
  #JOBDIRPREFIX=/vols/cms/rd1715/HiggsToInv/jobs_lighttree_170201_ICHEP
  JOBDIRPREFIX=/vols/cms/magnan/Hinvisible/RunIILT/jobs_lighttree_170209
  JOBDIR=$JOBDIRPREFIX/
  #OUTPUTPREFIX=/vols/cms/rd1715/HiggsToInv/output_lighttree_${PRODUCTION}_ICHEP
  #OUTPUTPREFIX=/vols/cms/rd1715/HiggsToInv/output_lighttree_170201_ICHEP
  OUTPUTPREFIX=/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170209

  OUTPUTDIR=$OUTPUTPREFIX/

  if [ "$SYST" != "central" ]
    then
    JOBDIR=$JOBDIRPREFIX/$SYST/
    OUTPUTDIR=$OUTPUTPREFIX/$SYST/
  fi

  if [ "$SYST" = "JESUP" ]
    then
    SYSTOPTIONS="--dojessyst=true --jesupordown=true --dojersyst=false"
  fi

  if [ "$SYST" = "JESDOWN" ]
    then
    SYSTOPTIONS="--dojessyst=true --jesupordown=false --dojersyst=false"
  fi

  if [ "$SYST" = "JERBETTER" ]
    then
    SYSTOPTIONS="--dojessyst=false --dojersyst=true --jerbetterorworse=true"
  fi

  if [ "$SYST" = "JERWORSE" ]
    then
    SYSTOPTIONS="--dojessyst=false --dojersyst=true --jerbetterorworse=false"
  fi

  if [ "$SYST" = "UESUP" ]
    then
    SYSTOPTIONS="--douessyst=true --uesupordown=true"
  fi

  if [ "$SYST" = "UESDOWN" ]
    then
    SYSTOPTIONS="--douessyst=true --uesupordown=false"
  fi

  if [ "$SYST" = "ELEEFFUP" ]
    then
    SYSTOPTIONS="--doidisoerr=true --doidisoerrmuore=false --doidisoerrupordown=true"
  fi

  if [ "$SYST" = "ELEEFFDOWN" ]
    then
    SYSTOPTIONS="--doidisoerr=true --doidisoerrmuore=false --doidisoerrupordown=false"
  fi

  if [ "$SYST" = "MUEFFUP" ]
    then
    SYSTOPTIONS="--doidisoerr=true --doidisoerrmuore=true --doidisoerrupordown=true"
  fi

  if [ "$SYST" = "MUEFFDOWN" ]
    then
    SYSTOPTIONS="--doidisoerr=true --doidisoerrmuore=true --doidisoerrupordown=false"
  fi


  echo "Config file: $CONFIG"
  mkdir -p $JOBDIR
  mkdir -p $OUTPUTDIR

  cp $CONFIG $OUTPUTDIR

  for QUEUEDIR in short medium #long
    do
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

    #PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/${PRODUSER}/${PRODUCTION}_MC
    PREFIX=root://gfe02.grid.hep.ph.ic.ac.uk:1095//store/user/${PRODUSER}/${PRODUCTION}_MC

    if [ "$PRODUCTION" = "Dec18" ]
      then
      PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/${PRODUSER}/${PRODUCTION}/MC
    fi

    for FILELIST in `ls filelists/$PRODUCTION/$QUEUEDIR/*_MC_Powheg-VBF*125.dat`
    #for FILELIST in `ls filelists/$PRODUCTION/$QUEUEDIR/*_MC_EWK*`
    #for FILELIST in `ls filelists/$PRODUCTION/$QUEUEDIR/*_MC_*NuNu*`
    #for FILELIST in `ls filelists/$PRODUCTION/$QUEUEDIR/*_MC_*`
      do
      echo "Processing files in "$FILELIST

      echo $FILELIST
      echo $FILELIST > tmp.txt


      MCOPTION="--mc=spring16_80X"
      #PREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/${PRODUSER}/${PRODUCTION}_MC
      PREFIX=root://gfe02.grid.hep.ph.ic.ac.uk:1095//store/user/${PRODUSER}/${PRODUCTION}_MC
      sed "s/filelists\/${PRODUCTION}\/$QUEUEDIR\/${PRODUCTION}_MC_//" tmp.txt > tmp2.txt
      INPUTPARAMS="filelists/$PRODUCTION/Params${PRODUCTION}.dat"

      echo $PREFIX


      #sed "s/filelists\/${PRODUCTION}\/$QUEUEDIR\/${PRODUCTION}_MC_//" tmp.txt > tmp2.txt

      JOB=MC_`sed "s/\.dat//" tmp2.txt`

      echo "JOB name = $JOB"

      #grep "JetsToLL" tmp.txt
      #if (( "$?" != 0 )); then
      # continue
      #fi

      JPTCUT=$JPTCUTVAL
      grep "Htoinv" tmp.txt
      if (( "$?" == 0 )); then
        JPTCUT=0
        MCOPTION="--mc=spring16_80X --donoskim=true"
      fi

      NEEDSSTREAM=0
      grep "EWKW" tmp.txt
      if (( "$?" == 0 )); then
        NEEDSSTREAM=1
      fi
      grep  "JetsToLNu" tmp.txt
      if (( "$?" == 0 )); then
        NEEDSSTREAM=1
      fi
      if (( "$NEEDSSTREAM" == 1 )); then
        for FLAVOUR in enu munu taunu
          do
          WJOB=$JOB"_"$FLAVOUR

          $JOBWRAPPER $JOBDIR $OUTPUTDIR "./bin/$MYEXEC --cfg=$CONFIG --prod="$PRODUCTION" $MCOPTION --filelist="$FILELIST" --input_prefix=$PREFIX --output_name=$JOB.root --output_folder=$OUTPUTDIR $SYSTOPTIONS --inputparams=$INPUTPARAMS --wstream=$FLAVOUR --jet1ptcut="$JPTCUT" --jet2ptcut="$JPTCUT" --jettype=$JETTYPE &> $JOBDIR/$WJOB.log" $JOBDIR/$WJOB.sh $GRIDSETUP
          if [ "$DOSUBMIT" = "1" ]; then 
            $JOBSUBMIT $JOBDIR/$WJOB.sh    
          else
            echo "$JOBSUBMIT $JOBDIR/$WJOB.sh"
          fi
        done

      else
        $JOBWRAPPER $JOBDIR $OUTPUTDIR "./bin/$MYEXEC --cfg=$CONFIG --prod="$PRODUCTION" $MCOPTION --filelist="$FILELIST" --input_prefix=$PREFIX --output_name=$JOB.root --output_folder=$OUTPUTDIR $SYSTOPTIONS --inputparams=$INPUTPARAMS --jet1ptcut="$JPTCUT" --jet2ptcut="$JPTCUT" --jettype=$JETTYPE &> $JOBDIR/$JOB.log" $JOBDIR/$JOB.sh  $GRIDSETUP
        if [ "$DOSUBMIT" = "1" ]; then 
          $JOBSUBMIT $JOBDIR/$JOB.sh
        else
          echo "$JOBSUBMIT $JOBDIR/$JOB.sh"
        fi
      fi
      rm tmp.txt tmp2.txt

    done


#Process bkg common with HiggsTautau (ONLY NEEDED IN Mar20 and Apr04)
    DOSHARED=false
    if [ "$PRODUCTION" = "Mar20" ]
      then
      DOSHARED=true
      SHAREDPREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/rlane/Feb20/MC_53X/
      FILELISTPREFIX=Feb20_MC_53X_
    else
      if [ "$PRODUCTION" = "Apr04" ]
        then
        DOSHARED=true
        SHAREDPREFIX=root://xrootd.grid.hep.ph.ic.ac.uk//store/user/${PRODUSER}/$PRODUCTION/MCtaushared/
        FILELISTPREFIX=Apr04_MCtaushared_
      fi
    fi

    for FILELIST in filelists/$PRODUCTION/$QUEUEDIR/${FILELISTPREFIX}*
      do
      if [ "$DOSHARED" = "true" ]
        then
        echo "Processing files in "$FILELIST

        echo $FILELIST > tmp.txt

        sed "s/filelists\/${PRODUCTION}\/$QUEUEDIR\/${FILELISTPREFIX}//" tmp.txt > tmp2.txt
        JOB=MC_`sed "s/\.dat//" tmp2.txt`

        echo "JOB name = $JOB"

        #grep "JetsToLL" tmp.txt
        #if (( "$?" != 0 )); then
        #  continue
        #fi

        grep  "JetsToLNu" tmp.txt
        if (( "$?" == 0 )); then
          for FLAVOUR in enu munu taunu
            do

            WJOB=$JOB"_"$FLAVOUR
		
            $JOBWRAPPER $JOBDIR $OUTPUTDIR "./bin/$MYEXEC --cfg=$CONFIG --prod="$PRODUCTION" --filelist="$FILELIST" --input_prefix=$SHAREDPREFIX --output_name=$WJOB.root --output_folder=$OUTPUTDIR $SYSTOPTIONS --inputparams=$INPUTPARAMS --wstream=$FLAVOUR --jet1ptcut="$JPTCUT" --jet2ptcut="$JPTCUT" --jettype=$JETTYPE &> $JOBDIR/$WJOB.log" $JOBDIR/$WJOB.sh $GRIDSETUP
            if [ "$DOSUBMIT" = "1" ]; then 
              $JOBSUBMIT $JOBDIR/$WJOB.sh
            else
              echo "$JOBSUBMIT $JOBDIR/$WJOB.sh"
            fi
          done
        else
          $JOBWRAPPER $JOBDIR $OUTPUTDIR "./bin/$MYEXEC --cfg=$CONFIG --prod="$PRODUCTION" --filelist="$FILELIST" --input_prefix=$SHAREDPREFIX --output_name=$JOB.root --output_folder=$OUTPUTDIR $SYSTOPTIONS --inputparams=$INPUTPARAMS --jet1ptcut="$JPTCUT" --jet2ptcut="$JPTCUT" --jettype=$JETTYPE &> $JOBDIR/$JOB.log" $JOBDIR/$JOB.sh $GRIDSETUP
          if [ "$DOSUBMIT" = "1" ]; then 
            $JOBSUBMIT $JOBDIR/$JOB.sh
          else
            echo "$JOBSUBMIT $JOBDIR/$JOB.sh"
          fi
        fi

        rm tmp.txt tmp2.txt
      fi

    done

  done

done
#if (( "$#" != "2" ))
#then
#echo ""
#    exit
#fi
