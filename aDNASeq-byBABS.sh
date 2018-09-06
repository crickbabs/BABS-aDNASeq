#!/bin/bash

POSITIONAL=()
OUTDIR=""

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	--outdir)
	    OUTDIR="$2"
	    shift # past argument
	    shift # past value
	    ;;
	--design)
	    DESIGN="--design $2"
	    if [[ "$2" = /* ]];
	    then
		DESIGN="--design $2"
	    else
		DESIGN="--design $PWD/$2"
	    fi
	    shift # past argument
	    shift # past value
	    ;;
	*)    # unknown option
	    POSITIONAL+=("$1") # save it in an array for later
	    shift # past argument
	    ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [[ $OUTDIR == "" ]];
then
    echo "ERROR: --outdir not set"
    exit
elif [ ! -d $OUTDIR ];
then
    echo "$OUTDIR does not exist, creating...."
    mkdir -p $OUTDIR
fi

# get wrapper script home 
BASH_HOME=$( dirname "${BASH_SOURCE[0]}" )
# check if absolute, if not make so
if [[ "$BASH_HOME" = /* ]];
then
    BASH_HOME=$BASH_HOME
else
    BASH_HOME=$PWD/$BASH_HOME
fi

cd $OUTDIR
nextflow run $BASH_HOME/aDNASeq-byBABS.nf $DESIGN ${POSITIONAL[@]}
cd -
