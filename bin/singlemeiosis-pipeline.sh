#! /bin/bash

#
# SINGLE MEIOSIS PIPELINE
#

# Authors: Delphine Charif <Delphine.Charif@versailles.inra.fr>, Joseph Tran <Joseph.Tran@versailles.inra.fr>

# Date: 2013-10-11

VERSION=dev

########################
# SECTION CONFIGURATION
#######################

# Include lib functions

PROD_PREFIX="/usr/local"
DEV_PREFIX="$(pwd)/.."
PREFIX=$DEV_PREFIX # TO BE CHANGED WHEN SWITCHING TO PROD

. $PROD_PREFIX/share/bash-common/lib/bash-common_lib.inc
. $PREFIX/share/singlemeiosis-pipeline/lib/singlemeiosis-pipeline_lib.inc

# CLI ARGUMENTS

ARGS=3
SeqFile1=$1
SeqFile2=$3
SampleName=$3

# SESSION VARIABLES

DATE=$(date '+%F_%Hh%Mm%Ss')
SESSION_ID=$(date '+%Y%M%d%H%M%S')
EXECUTED_COMMAND="$0 $*"

NAMESPACE="SINGLEMEIOSIS"
JOB_TAG=$NAMESPACE_$SampleName
SESSION_TAG=${JOB_TAG}_${USER}_${SESSION_ID}

WORKING_DIR=$(pwd)
SESSION_DIR=$SESSION_TAG
LOG_DIR=$SESSION_DIR/"log"
LOGFILE=${SESSION_TAG}.log

ERROR_TMP="/tmp/$(basename ${0%.*})_error_${SESSION_TAG}.log"


# CONFIG PARAMETERS LIST TO BE MOVED TO CONFIG FILE: cf ../share/singlemeiosis-pipeline/etc/singlemeiosis-pipeline_user.config

PIPELINE_SHARED=$PREFIX/share/$(basename ${0%.*})
PROD_PIPELINE_USER_CONFIG=$WORKING_DIR/$(basename ${0%.*})_user.config
DEV_PIPELINE_USER_CONFIG=$PIPELINE_SHARED/etc/$(basename ${0%.*})_user.config
PIPELINE_USER_CONFIG=$DEV_PIPELINE_USER_CONFIG # TO BE CHANGED WHEN SWITCHING TO PROD

# GET GENOME ALIASES LIST
GENOME_ALIASES_LIST=($(get_genome_aliases_list_wo_snpeff $PIPELINE_USER_CONFIG 2>$ERROR_TMP))
# TODO: print a warning if error

#===============
# USAGE MESSAGE
#===============
[[ $# -ne "$ARGS" ]] && { printf %s "\
Program: $(basename $0)
Version: $VERSION
Contact: IJPB Bioinformatics Dev Team

Usage: $(basename $0) SeqFile1 SeqFile2 SampleName

Arguments: SeqFile1 	Forward read sequences file (Illumina fastq file)
           SeqFile2 	Reverse read sequences file (Illumina fastq file)
           SampleName  	Prefix to use for the analysis, i.e. prefix can be the sample name or anything else suitable

Genome aliases list to be used in the user config file (get a copy here: $DEV_PIPELINE_USER_CONFIG):
$(for ga in "${GENOME_ALIASES_LIST[@]}"; do
	echo "- $ga"
done)

Notes: 1. this pipeline version is actually able to perform reads mapping using bwa (version $(get_tool_version 'bwa')) and per base coverage computation using samtools (version $(get_tool_version 'samtools'))

";
exit 1; }


########################
# PIPELINE STEPS
#######################

## CREATE SESSION DIR
## LOAD CONFIG
## SET GENOME PATH AND INDEXES
## MAPPING
## FILTERING
## STATISTICS
## ANALYSIS
## CLEANING

























