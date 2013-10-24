#! /bin/bash

#
# SINGLE MEIOSIS PIPELINE
#

# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>, Delphine Charif <Delphine.Charif@versailles.inra.fr>

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

MandatoryArgs=2
ConfigFile=$1
SessionName=$2
EmailAdress=$3

# SESSION VARIABLES

DATE=$(date '+%F_%Hh%Mm%Ss')
SESSION_ID=$(date '+%Y%M%d%H%M%S')
EXECUTED_COMMAND="$0 $*"

NAMESPACE="SINGLEMEIOSIS"
SESSION_TAG=${SessionName}_${USER}_${SESSION_ID}

WORKING_DIR=$(pwd)
SESSION_DIR=$SessionName
LOG_DIR=$SESSION_DIR/"log"
LOGFILE=${SESSION_TAG}.log

ERROR_TMP="/tmp/$(basename ${0%.*})_error_${SESSION_TAG}.log"


# CONFIG PARAMETERS LIST TO BE MOVED TO CONFIG FILE: cf ../share/singlemeiosis-pipeline/etc/singlemeiosis-pipeline_user.config

PIPELINE_SHARED=$PREFIX/share/$(basename ${0%.*})
DEV_PIPELINE_USER_CONFIG=$PIPELINE_SHARED/etc/$(basename ${0%.*})_user.config
PIPELINE_USER_CONFIG=$DEV_PIPELINE_USER_CONFIG
BACKUPED_CONFIG_FILE=$SESSION_DIR/$(basename $ConfigFile)

# GET GENOME ALIASES LIST
GENOME_ALIASES_LIST=($(get_genome_aliases_list_wo_snpeff $PIPELINE_USER_CONFIG 2>$ERROR_TMP))
rtrn=$?
genome_alias_list_failed_msg="[Genome alias list] Failed An error occured while retrieving the list of genome alias. See $ERROR_TMP file for more details."
if [[ "$rtrn" -ne 0 ]]; then
	echo -e "$(date '+%Y_%m_%d %T') $genome_alias_list_failed_msg" 2>&1
	exit $rtrn
fi

#===============
# USAGE MESSAGE
#===============
[[ $# -ne "$MandatoryArgs" ||  $# -ne $(expr $MandatoryArgs + 1) ]] && { printf %s "\
Program: $(basename $0)
Version: $VERSION
Contact: IJPB Bioinformatics Dev Team

Usage: $(basename $0) ConfigFile SessionName [EmailAdress]

Arguments: ConfigFile     The user configuration file listing the tetrad analysis parameters and data samples paths
                          You can get a copy there: $PIPELINE_USER_CONFIG
                          After copying this file in your working directory, edit this file to map your needs.
           SessionName    Prefix to use for the tetrad analysis, i.e. prefix can be the tetrad name or anything else suitable without spaces
           [EmailAdress]  An optional but valid email address to send pipeline job/error status notifications

Available genome aliases to be used in config file are:
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

#===================
# SESSION DIRECTORY 
#===================

echo "$(date '+%Y-%m-%d %T') [$(basename $0)] Start running $NAMESPACE pipeline (version: $VERSION)." | tee $ERROR_TMP 2>&1
echo "$(date '+%Y-%m-%d %T') [$(basename $0)] Executed command: $0 $*" | tee -a $ERROR_TMP 2>&1

#
# Create a directory named with SESSION_DIR value, to save all outputs 
#
echo "$(date '+%Y-%m-%d %T') [Session directory] Creating $SESSION_DIR directory ..." | tee -a $ERROR_TMP 2>&1
if [[ -d $SESSION_DIR ]]; then
    echo "$(date '+%Y-%m-%d %T') [Session directory] OK $SESSION_DIR directory already exists. Will output all session files in this directory." | tee -a $ERROR_TMP 2>&1
else
    mkdir $SESSION_DIR 2>>$ERROR_TMP
	rtrn=$?
	job_dir_failed_msg="[Session directory] Failed Session directory, $SESSION_DIR, was not created."
	exit_on_error "$ERROR_TMP" "$job_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" $SESSION_TAG $EmailAdress
	echo "$(date '+%Y-%m-%d %T') [Job directory] OK $SESSION_DIR directory was created successfully. Will output all session files in this directory." | tee -a $ERROR_TMP 2>&1
fi

# Create log directory
echo "$(date '+%Y-%m-%d %T') [Log directory] Creating $LOG_DIR directory ..." | tee -a $ERROR_TMP 2>&1
if [[ -d $LOG_DIR ]]; then
    [[ -s $ERROR_TMP ]] && cat $ERROR_TMP > $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y-%m-%d %T') [Log directory] OK $LOG_DIR directory already exists. Will write log files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    mkdir $LOG_DIR 2>>$ERROR_TMP
	rtrn=$?
	log_dir_failed_msg="[Log directory] Failed Log directory, $LOG_DIR, was not created."
	exit_on_error "$ERROR_TMP" "$log_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" $SESSION_TAG $EmailAdress 
	[[ -s $ERROR_TMP ]] && cat $ERROR_TMP > $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y-%m-%d %T') [Log directory] OK $LOG_DIR directory was created sucessfully. Will write log files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1	
fi

#=============
# LOAD CONFIG 
#=============

#
# Test for absence of user config file
# if exists ok continue else display error message and exit
#
echo "$(date '+%Y-%m-%d %T') [Check config: session user config file] Checking if user config file $ConfigFile exists ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ -s $ConfigFile ]]; then
    echo "$(date '+%Y-%m-%d %T') [Check config: session user config file] OK User config file, $ConfigFile, exists and is not empty." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
	user_config_failed_msg="[Check config: session user config file] Failed User config file, $ConfigFile, does not exist or is empty."
	echo -e "$user_config_failed_msg" 2>&1 >$ERROR_TMP
    echo "$(date '+%Y-%m-%d %T') [Check config: session user config file] Error: " | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit_on_error "$ERROR_TMP" "$user_config_failed_msg" 3 "$LOG_DIR/$LOGFILE" $SESSION_TAG $EmailAdress
fi

# 1. Backup session user config file in session dir if not exist
echo "$(date '+%Y-%m-%d %T') [Check config: session user config file] Backuping session user config file into session directory ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
cp $ConfigFile $SESSION_DIR/. 2>$ERROR_TMP
rtrn=$?
cp_user_config_failed_msg="[Check config: session user config file] Failed backuping session user config file into session directory."
exit_on_error "$ERROR_TMP" "$cp_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" $SESSION_TAG $EmailAdress
echo "$(date '+%Y-%m-%d %T') [Check config: session user config file] Will use backuped session user config file: $BACKUPED_CONFIG_FILE" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Load config parameters from backuped session user config file
echo "$(date '+%Y-%m-%d %T') [Check config: session user config file] Loading session user config parameters from $BACKUPED_CONFIG_FILE file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
load_user_config_failed_msg="[Check config: session user config file] Failed loading session user config parameters from $BACKUPED_CONFIG_FILE file."
for cfg in $(get_config_sections $BACKUPED_CONFIG_FILE 2>$ERROR_TMP;); do
	rtrn=$?	
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" $SESSION_TAG $EmailAdress
    echo -e "--- Config section [${cfg}] ---"
    unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }' 2>$ERROR_TMP) $(toupper ${NAMESPACE}_${cfg}_) 2>>$ERROR_TMP
	rtrn=$?
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" $SESSION_TAG $EmailAdress
    set_config_params $BACKUPED_CONFIG_FILE ${cfg} ${NAMESPACE} 2>$ERROR_TMP
    rtrn=$?
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" $SESSION_TAG $EmailAdress
    for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
		echo -e "$params"
    done
	rtrn=$?
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" $SESSION_TAG $EmailAdress
done
echo "$(date '+%Y-%m-%d %T') [Check config: session user config file] OK Session user config file, $BACKUPED_CONFIG_FILE, was loaded successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1


























