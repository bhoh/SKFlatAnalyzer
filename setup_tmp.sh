export SKFlat_WD=`pwd`
export SKFlat_LIB_PATH=$SKFlat_WD/lib/
mkdir -p $SKFlat_LIB_PATH
mkdir -p $SKFlat_WD/tar

export SKFlatV="Run2Legacy_v1"
mkdir -p $SKFlat_WD/data/$SKFlatV
export DATA_DIR=$SKFlat_WD/data/$SKFlatV

#### USER INFO ####
export SKFlatLogEmail=''
export SKFlatLogWeb='' # leave if blank if you don't have webpage
export SKFlatLogWebDir='' # leave if blank if you don't have webpage

#### use cvmfs for root ####
export CMS_PATH=/cvmfs/cms.cern.ch
source $CMS_PATH/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
export cmsswrel='cmssw-patch/CMSSW_10_4_0_patch1'
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/src
echo "@@@@ SCRAM_ARCH = "$SCRAM_ARCH
echo "@@@@ cmsswrel = "$cmsswrel
echo "@@@@ scram..."
eval `scramv1 runtime -sh`
cd -
source /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/external/$SCRAM_ARCH/bin/thisroot.sh

if [[ $HOSTNAME == *"ui"*".sdfarm.kr"* ]]; then

  echo "@@@@ Working on KISTI"
  export SKFlatRunlogDir="/cms/ldap_home/$USER/SKFlatRunlog/"
  export SKFlatOutputDir="/cms/ldap_home/$USER/SKFlatOutput/"

elif [[ $HOSTNAME == *"cms.snu.ac.kr"* ]]; then

  echo "@@@@ Working on tamsa1"
  export SKFlatRunlogDir="/data7/Users/$USER/SKFlatRunlog/"
  export SKFlatOutputDir="/data7/Users/$USER/SKFlatOutput/"

elif [[ $HOSTNAME == *"tamsa2"* ]]; then

  echo "@@@@ Working on tamsa2"
  export SKFlatRunlogDir="/data7/Users/$USER/SKFlatRunlog/"
  export SKFlatOutputDir="/data7/Users/$USER/SKFlatOutput/"

elif [[ $HOSTNAME == *"knu"* ]]; then

  echo "@@@@ Working on KNU"
  export SKFlatRunlogDir="/u/user/$USER/scratch/SKFlatRunlog/"
  export SKFlatOutputDir="/u/user/$USER/scratch/SKFlatOutput/"

fi

alias skout="cd $SKFlatOutputDir/$SKFlatV/"

export MYBIN=$SKFlat_WD/bin/
export PYTHONDIR=$SKFlat_WD/python/
export PATH=${MYBIN}:${PYTHONDIR}:${PATH}

export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SKFlat_WD/DataFormats/include/:$SKFlat_WD/AnalyzerTools/include/:$SKFlat_WD/Analyzers/include/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKFlat_LIB_PATH

source $SKFlat_WD/bin/BashColorSets.sh

## Todo list ##
python python/PrintToDoLists.py
source $SKFlat_WD/tmp/ToDoLists.sh
rm $SKFlat_WD/tmp/ToDoLists.sh

## Log Dir ##
echo "* Your Log Directory Usage"
du -sh $SKFlatRunlogDir
echo "-----------------------------------------------------------------"
