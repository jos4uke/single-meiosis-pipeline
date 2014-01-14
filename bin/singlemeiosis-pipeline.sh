#! /bin/bash

#
# SINGLE MEIOSIS PIPELINE
#

# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>, Delphine Charif <Delphine.Charif@versailles.inra.fr>

# This script provides a pipeline for processing single meiosis data

# This software is governed by the CeCILL license, Version 2.0 (the "License"), under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license, Version 2.0 (the "License"), as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt". 

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license, Version 2.0 (the "License"), and that you accept its terms.

# Date: 2013-10-11

declare -r VERSION="dev"


########################
# SECTION CONFIGURATION
#######################

### SESSION VARIABLES ###

NAMESPACE="SINGLEMEIOSIS"

WORKING_DIR=$(pwd)
DATE=$(date '+%F_%Hh%Mm%Ss')
SESSION_ID=$(date '+%Y%M%d%H%M%S')
EXECUTED_COMMAND="$0 $*"
SESSION_TAG=${NAMESPACE}_${USER}_${SESSION_ID}

LOG_DIR="log"
DEBUGFILE=${SESSION_TAG}.log
ERROR_TMP="/tmp/$(basename ${0%.*})_error_${SESSION_TAG}.log"

[[ $VERSION -eq "dev" ]] && PROG_PATH=$(realpath $(dirname $0));PIPELINE_USER_CONFIG=${PROG_PATH}/../share/singlemeiosis-pipeline/etc/singlemeiosis-pipeline_user.config || PIPELINE_USER_CONFIG=/usr/local/share/singlemeiosis-pipeline/etc/singlemeiosis-pipeline_user.config

PIDS_ARR=()
WAITALL_TIMEOUT=259200
WAITALL_INTERVAL=60
WAITALL_DELAY=60

### LOGGING CONFIGURATION ###

# load log4sh (disabling properties file warning) and clear the default
# configuration
LOG4SH_CONFIGURATION='none' . /usr/local/share/log4sh/build/log4sh 2>/dev/null
[[ $? != 0 ]] && $(echo "Error loading log4sh lib" >&2; exit 1)
log4sh_resetConfiguration

# set the global logging level
logger_setLevel DEBUG

# add and configure a FileAppender that outputs to STDERR
logger_addAppender stderr
appender_setType stderr FileAppender
appender_file_setFile stderr STDERR
appender_setLevel stderr FATAL
appender_setLayout stderr PatternLayout
appender_setPattern stderr '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m' 
appender_activateOptions stderr
appender_exists stderr && logger_debug "Standard error appender is enabled." || logger_warn "Standard error appender was not enabled. Maybe a log4sh error occured."

# add and configure console appender that outputs to standard output
logger_addAppender console
appender_setType console ConsoleAppender
appender_setLevel console INFO 
appender_setLayout console PatternLayout
appender_setPattern console '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m' 
appender_activateOptions console
appender_exists console && logger_debug "Console appender is enabled." || logger_warn "Console appender was not enabled. Maybe a log4sh error occured."


### LOAD LIB ###

# bash-common lib
[[ $VERSION -eq "dev" ]] && LIB_PATH=$(realpath $(dirname $0))/../../bash-common/share/bash-common/lib/bash-common_lib.inc || LIB_PATH=/usr/local/share/bash-common/lib/bash-common_lib.inc

. $LIB_PATH
if [[ $? -ne 0 ]]; then
	logger_fatal "Error loading bash common lib: $LIB_PATH"  
	exit 1
fi

# singlemeiosis lib
[[ $VERSION -eq "dev" ]] && LIB_PATH=$(realpath $(dirname $0))/../share/singlemeiosis-pipeline/lib/singlemeiosis-pipeline_lib.inc || LIB_PATH=/usr/local/share/singlemeiosis-pipeline/lib/singlemeiosis-pipeline_lib.inc

. $LIB_PATH
if [[ $? -ne 0 ]]; then 
	logger_fatal "Error loading singlemeiosis lib: $LIB_PATH"
	exit 1
fi


### GET GENOME ALIASES LIST ###
GENOME_ALIASES_LIST=($(get_genome_aliases_list_w_bwa $PIPELINE_USER_CONFIG 2>$ERROR_TMP))
rtrn=$?
genome_alias_list_failed_msg="[Genome alias list] Failed An error occured while retrieving the list of genome alias. See $ERROR_TMP file for more details."
if [[ "$rtrn" -ne 0 ]]; then
	logger_fatal "$genome_alias_list_failed_msg"
	exit $rtrn
fi



### USAGE ###
Usage()
{
printf %s "\
Program: $(basename $0)
Version: $VERSION

Copyright 2013 Joseph Tran <Joseph.Tran@versailles.inra.fr> & Delphine Charif <Delphine.Charif@versailles.inra.fr>

Licensed under the CeCILL License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

Usage: $(basename $0) -c|--configfile CONFIG_FILE -o|--out_dir OUTPUT_DIR [-d|--debug] [-e|--email_address VALID_EMAIL_ADDR]

Options:
-c|--config_file CONFIG_FILE            The user configuration file listing the data samples paths and tetrad analysis parameters.
                                        You can get a copy there: $PIPELINE_USER_CONFIG.
                                        See 'Configuration file format' section below.
-o|--out_dir OUTPUT_DIR                 The output directory.
-d|--debug                              Enable debugging mode in the console. 
-e|--email_address VALID_EMAIL_ADDR     An optional but valid email address to send pipeline job/error status notifications
-h|--help                               Displays this message.

Configuration file format:
TODO: describe the config section format
Note that blank lines are ignored.

[config_section_1]
variable1=value1
variable2=value2
# this is a comment line: put a description here
[config_section_1]
variable1=value1
variable2=value2

Available genome aliases to be used in genome_alias section in config file are:
$(for ga in "${GENOME_ALIASES_LIST[@]}"; do
	echo "- $ga"
done)

example:

[genome_alias]
# parent 1
papa=Ath-Col0-Vers-WG
# parent 2
mama=Ath-Ler0-Vers-WG
# parents 1 & 2
papamama=Ath-Coller0_Vers-WG

Notes: 1. this pipeline version is actually able to perform reads mapping using bwa (version $(get_tool_version 'bwa')) and per base coverage computation using samtools (version $(get_tool_version 'samtools'))

"
}

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; 
CONFIGURE_OPTS=`getopt -o hc:o:e:d --long help,config_file:,out_dir:,debug,email_address: \
	-n 'singlemeiosis-pipeline.sh' -- "$@"`

if [[ $? != 0 ]] ; then Usage >&2 ; exit 1 ; fi

# Note the quotes around `$CONFIGURE_OPTS'
eval set -- "$CONFIGURE_OPTS"

while true; do
	case "$1" in
		-h | --help ) Usage >&2; exit 1;;
		-c | --config_file ) CONFIGFILE="$2"; shift 2 ;;
		-o | --out_dir ) OUTPUT_DIR="$2"; shift 2 ;;
		-d | --debug ) 
					appender_setLevel console DEBUG;
					appender_activateOptions console;
					shift 1 ;;
		-e | --email_address ) EMAIL="$2"; shift 2 ;;
		-- ) shift; break ;;
		* ) break ;;
	esac
done


### VALIDATION ###
if [[ ! -s $CONFIGFILE ]]; then 
	logger_fatal "Config file, $CONFIGFILE, does not exist or is empty. See Usage with --help option."; 
	exit 1;
fi
# TODO: email validator

########################
# PIPELINE STEPS
#######################

## CREATE OUTPUT DIR
## LOAD CONFIG
## SET GENOME PATH AND INDEXES
## MAPPING
## FILTERING
## STATISTICS
## ANALYSIS
## CLEANING

#===================
# OUTPUT DIRECTORY 
#===================

echo "Start running $NAMESPACE pipeline (version: $VERSION)." | tee $ERROR_TMP 2>&1 | logger_info
echo "Executed command: $0 $*" | tee -a $ERROR_TMP 2>&1 | logger_info

#
# Create a directory named with OUTPUT_DIR value, to save all outputs 
#
echo "Creating $OUTPUT_DIR directory ..." | tee -a $ERROR_TMP 2>&1 | logger_info
if [[ -d $OUTPUT_DIR ]]; then
    echo "OK $OUTPUT_DIR directory already exists. Will output all output files in this directory." | tee -a $ERROR_TMP 2>&1  | logger_info
else
    mkdir $OUTPUT_DIR 2>>$ERROR_TMP
	rtrn=$?
	out_dir_failed_msg="[Output directory] Failed. Output directory, $OUTPUT_DIR, was not created."
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$out_dir_failed_msg"
	exit_on_error "$ERROR_TMP" "$out_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
	echo "$(date '+%Y-%m-%d %T') [Output directory] OK $OUTPUT_DIR directory was created successfully. Will output all output files in this directory." | tee -a $ERROR_TMP 2>&1 | logger_info
fi

# Create log directory
echo "Creating $LOG_DIR directory ..." | tee -a $ERROR_TMP 2>&1 | logger_info
if [[ -d $OUTPUT_DIR/$LOG_DIR ]]; then
    echo "OK $OUTPUT_DIR/$LOG_DIR directory already exists. Will write log files in this directory." | tee -a $ERROR_TMP 2>&1 | logger_info
else
    mkdir $OUTPUT_DIR/$LOG_DIR 2>>$ERROR_TMP
	rtrn=$?
	log_dir_failed_msg="[Log directory] Failed Log directory, $OUTPUT_DIR/$LOG_DIR, was not created."
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$log_dir_failed_msg"
	exit_on_error "$ERROR_TMP" "$log_dir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL 
	echo "$(date '+%Y-%m-%d %T') [Log directory] OK $OUTPUT_DIR/$LOG_DIR directory was created sucessfully. Will write log files in this directory." | tee -a $ERROR_TMP 2>&1 | logger_info	
fi

#==================================
# Enable the pipeline debug logger
#==================================

logger_addAppender debuggerF 
appender_setType debuggerF FileAppender
appender_file_setFile debuggerF $OUTPUT_DIR/$LOG_DIR/$DEBUGFILE
appender_setLevel debuggerF DEBUG 
appender_setLayout debuggerF PatternLayout
appender_setPattern debuggerF '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m' 
appender_activateOptions debuggerF
appender_exists debuggerF && cat $ERROR_TMP | logger_info
appender_exists debuggerF && logger_info "Debugging infos will be output to $OUTPUT_DIR/$LOG_DIR/$DEBUGFILE file." || logger_warn "The debugger file appender was not enabled. Maybe a log4sh error occured."


#=============
# LOAD CONFIG 
#=============

# set backup config file variable
BACKUPED_CONFIG_FILE=$OUTPUT_DIR/$(basename $CONFIGFILE)

# 1. Backup session user config file in session dir if not exist
logger_info "[Check config: session user config file] Backuping session user config file into session directory ..."
cp $CONFIGFILE $OUTPUT_DIR/. 2>$ERROR_TMP
rtrn=$?
cp_user_config_failed_msg="[Check config: session user config file] Failed backuping session user config file into session directory."
[[ "$rtrn" -ne 0 ]] && logger_fatal "$cp_user_config_failed_msg"
exit_on_error "$ERROR_TMP" "$cp_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
logger_info "[Check config: session user config file] Will use backuped session user config file: $BACKUPED_CONFIG_FILE" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Load config parameters from backuped session user config file
logger_info "[Check config: session user config file] Loading session user config parameters from $BACKUPED_CONFIG_FILE file ..."
load_user_config_failed_msg="[Check config: session user config file] Failed loading session user config parameters from $BACKUPED_CONFIG_FILE file."
for cfg in $(get_config_sections $BACKUPED_CONFIG_FILE 2>$ERROR_TMP;); do
	rtrn=$?	
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
    logger_debug "--- Config section [${cfg}] ---"
    unset $(set |awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
		  pattern = "\^" prefix "_" cfg "_" 
       }
       $0~pattern { print $1 }' 2>$ERROR_TMP ) 2>>$ERROR_TMP
	rtrn=$?
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
    CONFIG_PARAMS=$(format_config_params $BACKUPED_CONFIG_FILE ${cfg} ${NAMESPACE} 2>$ERROR_TMP)
	eval "${CONFIG_PARAMS}"
    rtrn=$?
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
    for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
		logger_debug "$params"
    done
	rtrn=$?
	[[ "$rtrn" -ne 0 ]] && logger_fatal "$load_user_config_failed_msg"
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
done
logger_info "[Check config: session user config file] OK Session user config file, $BACKUPED_CONFIG_FILE, was loaded successfully."

#==========================================
# GENOMES AND INDEX PATH
#==========================================

## SET GENOME FASTA FILES PATH

logger_info "[Genome sequences and index paths] set path variables ..."
### set genome fasta paths for the papa and mama genomes
declare -r genome_base_path=$(toupper ${NAMESPACE}_paths)_GENOMES_BASE_PATH

if [[ -z ${!genome_base_path} && ! -d ${!genome_base_path} ]]; then 
	logger_fatal "An error occured while setting genome base path variable."
	exit 1
fi
logger_debug "[Genome base path] ${genome_base_path}=${!genome_base_path}"

#### papa
declare -r ga_papa=$(toupper ${NAMESPACE}_genome_alias )_papa
if [[ -z ${!ga_papa} ]]; then 
	logger_fatal "An error occured while setting genome alias variable for papa genome."
	exit 1
fi
logger_debug "[Genome alias] ${ga_papa}=${!ga_papa}"

eval "$(toupper ${NAMESPACE}_paths)_papa_fasta=${!genome_base_path}/${!ga_papa}/$(ls ${!genome_base_path}/${!ga_papa} | grep -e "${!ga_papa}\.m*fas*$")"
declare -r ga_papa_fasta=$(toupper ${NAMESPACE}_paths)_papa_fasta
if [[ ! -s ${!ga_papa_fasta} ]]; then 
	logger_fatal "An error occured while setting genome alias sequence path for papa genome."
	exit 1
fi
logger_info "[Genome alias sequence path] ${ga_papa_fasta}=${!ga_papa_fasta}"

# call papa directly
#eval echo -e \$"$(toupper ${NAMESPACE}_paths)_papa_fasta"

##### mama
declare -r ga_mama=$(toupper ${NAMESPACE}_genome_alias )_mama
if [[ -z ${!ga_mama} ]]; then 
	logger_fatal "An error occured while setting genome alias variable for mama genome."
	exit 1
fi
logger_debug "[Genome alias] ${ga_mama}=${!ga_mama}"

eval "$(toupper ${NAMESPACE}_paths)_mama_fasta=${!genome_base_path}/${!ga_mama}/$(ls ${!genome_base_path}/${!ga_mama} | grep -e "${!ga_mama}\.m*fas*$")"
declare -r ga_mama_fasta=$(toupper ${NAMESPACE}_paths)_mama_fasta
if [[ ! -s ${!ga_mama_fasta} ]]; then 
	logger_fatal "An error occured while setting genome alias sequence path for mama genome."
	exit 1
fi
logger_info "[Genome alias sequence path] ${ga_mama_fasta}=${!ga_mama_fasta}"

# call mama directly
#eval echo -e \$"$(toupper ${NAMESPACE}_paths)_mama_fasta"

## SET GENOME BWA INDEX PATH RELATIVE TO CURRENT VERSION/TOOL

#### set current tool version index for the papamama genome
declare -r genome_index_path=$(toupper ${NAMESPACE}_paths)_INDEXES_BASE_PATH
if [[ -z ${!genome_index_path} ]]; then 
	logger_fatal "An error occured while setting genome indexes path variable."
	exit 1
fi
logger_debug "[Genome index path] ${genome_index_path}=${!genome_index_path}"

declare -r genome_bwa_path=$(toupper ${NAMESPACE}_paths)_BWA_INDEXES
if [[ -z ${!genome_bwa_path} ]]; then 
	logger_fatal "An error occured while setting genome bwa indexes path variable."
	exit 1
fi
logger_debug "[Genome index path] ${genome_bwa_path}=${!genome_bwa_path}"

declare -r ga_papamama=$(toupper ${NAMESPACE}_genome_alias)_papamama
if [[ -z ${!ga_papamama} ]]; then 
	logger_fatal "An error occured while setting genome alias variable for the papamama genome."
	exit 1
fi
logger_debug "[Genome alias] ${ga_papamama}=${!ga_papamama}"

ext="fas" # TODO: remove fas extension from index files, keep only genome alias
eval "$(toupper ${NAMESPACE}_paths)_papamama_bwa_index=${!genome_index_path}/${!genome_bwa_path}/$(get_tool_version bwa)/${!ga_papamama}/${!ga_papamama}.$ext"
declare -r papamama_bwa_index_path=$(toupper ${NAMESPACE}_paths)_papamama_bwa_index
if [[ -z ${!papamama_bwa_index_path} ]]; then 
	logger_fatal "An error occured while setting genome bwa index path variable for the papamama genome."
	exit 1
fi
IDX_FILES=($(ls ${!papamama_bwa_index_path}*))
if [[ ${#IDX_FILES[@]} -le 0 ]]; then
	logger_fatal "An error occured while checking genome bwa index files for the papamama genome."
	exit 1
fi
logger_info "[Genome index path] ${papamama_bwa_index_path}=${!papamama_bwa_index_path}"

# call directly
#eval echo -e \$"$(toupper ${NAMESPACE}_paths)_papamama_bwa_index"
# test
#eval ls -lh \$"$(toupper ${NAMESPACE}_paths)_papamama_bwa_index*"


#==========================================
# TETRAD SAMPLES
#==========================================

# F1 & M1-4, etc.
tscs="tetrad_samples"
TETRAD_SAMPLES=($(set |awk -F= -vcfg="${tscs}" -vpfx="${NAMESPACE}" -vsfx="sample_name_alias" 'BEGIN { 
          cfg = toupper(cfg);
          pfx = toupper(pfx);
			pattern = pfx "_" cfg "_\.\*_" sfx;
       }
       $0~pattern { print $1 }' 2>$ERROR_TMP))
logger_info "[Tetrad samples] ${#TETRAD_SAMPLES[@]} samples in $tscs section loaded from config file."

logger_info "[Tetrad samples] check for sequence files for each tetrad sample"
SAMPLES_STACK=()
SKIPPED=()
for s in "${TETRAD_SAMPLES[@]}"; do
	logger_debug "$s=${!s}"
	sid=$(echo $s | awk -F"_" '{print $4}')
	logger_debug "sample id: $sid"
	seqFR=($(set | grep -e "$(toupper ${NAMESPACE}_${tscs})_$sid_.*_seqfile_R" | cut -d\= -f1))
	logger_debug "${seqFR[0]}=${!seqFR[0]}"
	logger_debug "${seqFR[1]}=${!seqFR[1]}"
	if [[ -s  ${!seqFR[0]} && -s ${!seqFR[1]} ]]; then
		logger_info "The given pair of fastq files does exist for the current sample ${!s}. Add the sample to the stack."
		SAMPLES_STACK=("${SAMPLES_STACK[@]}" "$s")
	else
		[[ ! -s "${!seqFR[0]}" ]] && logger_warn "${seqFR[0]}=${!seqFR[0]} file does not exist or is empty."
		[[ ! -s "${!seqFR[1]}" ]] && logger_warn "${seqFR[1]}=${!seqFR[1]} file does not exist or is empty."
		logger_warn "No pair of fastq files does exist for the current sample ${!s}. Skip this sample as it will not be added to the stack."
		SKIPPED=("${SKIPPED[@]}" "$s")
	fi
done

if [[ "${#SAMPLES_STACK}" -eq 0 ]]; then
	logger_warn "The samples stack is empty."
	logger_warn "All the samples were skipped."
	logger_debug "Skipped samples: ${SKIPPED[@]}"
	logger_info "Exit the pipeline."
	exit 0
else
	logger_info "The samples stack contains ${#SAMPLES_STACK[@]} samples to be processed."
	logger_info "Samples in the stack:"
	for s in "${SAMPLES_STACK[@]}"; do
		logger_info "$s=${!s}"
	done
	logger_info "Samples skipped:"
	for sk in "${SKIPPED[@]}"; do
		logger_info "$sk=${!sk}"
	done	
fi

#==========================================
# MAPPING
#==========================================
logger_info "[Mapping] Run bwa on stack samples."
mapping_cmd=

# bwa aln
mapping_cmd="bwa aln"
for s in "${SAMPLES_STACK[@]}"; do
	logger_info "[Mapping] Run bwa aln on stack samples."
	logger_info "[Mapping] Current sample: ${!s}"
	# create the sample output dir $OUTPUT_DIR/${!s}
	if [[ ! -d $OUTPUT_DIR/${!s} ]]; then	
		mkdir $OUTPUT_DIR/${!s} 2>$ERROR_TMP; rtrn=$?
		if [[ $rtrn -ne 0 ]]; then 
			mkdir_err_msg="[Mapping] An error occured while creating $OUTPUT_DIR/${!s} directory."
			logger_fatal "$mkdir_err_msg"
			exit_on_error "$ERROR_TMP" "$mkdir_err_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
	fi
	logger_info "[Mapping] bwa will output files in $OUTPUT_DIR/${!s} directory." 

	# fastq files
	sid=$(echo $s | awk -F"_" '{print $4}')
	logger_debug "sample id: $sid"
	seqFR=($(set | grep -e "$(toupper ${NAMESPACE}_${tscs})_$sid_.*_seqfile_R" | cut -d\= -f1))

	for seqF in "${seqFR[@]}"; do
		logger_debug "[Mapping] Current fastq file: ${!seqF}"
		# error logging
		CURRENT_MAPPING_ERROR=$OUTPUT_DIR/${!s}/${!seqF}_mapping_err.log

		# build cli options
		bwa_aln_cli_options=($(buildCommandLineOptions "$mapping_cmd" "$NAMESPACE" 2>$CURRENT_MAPPING_ERROR))
		rtrn=$?
		cli_options_failed_msg="[Mapping] An error occured while building the $mapping_cmd command line options for current sample ${!s} and current fastq file ${!seqF}."
		exit_on_error "$CURRENT_MAPPING_ERROR" "$cli_options_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
		opts="${bwa_aln_cli_options[@]}"
		logger_debug "[Mapping] $mapping_cmd options: $opts"

		# build cli
		bwa_aln_cli="$mapping_cmd $opts ${!papamama_bwa_index_path} ${!seqF} >$OUTPUT_DIR/${!s}/${!seqF%.*}.sai 2>$CURRENT_MAPPING_ERROR &"

		# run the cli
		logger_debug "[Mapping] $bwa_aln_cli"
		eval "$bwa_aln_cli" 2>$ERROR_TMP
		pid=$!
		rtrn=$?
		eval_failed_msg="[Mapping] An error occured while eval $mapping_cmd cli."
		exit_on_error "$ERROR_TMP" "$eval_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
		logger_debug "[Mapping] $mapping_cmd pid: $pid"

		# add pid to array
		PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
	done
done

# wait until all bwa aln processes to finish then proceed to next step
# reinit pid array
pid_list_failed_msg="[Mapping] Failed getting process status for process $p."	
for p in "${PIDS_ARR[@]}"; do
	logger_debug "$(ps aux | grep $USER | gawk -v pid=$p '$2 ~ pid {print $0}' 2>${ERROR_TMP})"
	rtrn=$?
	exit_on_error "$ERROR_TMP" "$pid_list_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
done
logger_info "[Mapping] Wait for all $mapping_cmd processes to finish before proceed to next step."
waitalluntiltimeout "${PIDS_ARR[@]}" 2>/dev/null
logger_info "[Mapping] All $mapping_cmd processes finished. Will proceed to next step: bwa sampe."
PIDS_ARR=()

# bwa sampe
mapping_cmd="bwa sampe"
	logger_info "[Mapping] Run bwa aln on stack samples."
for s in "${SAMPLES_STACK[@]}"; do
	logger_info "[Mapping] Current sample: ${!s}"
	
	# fastq files
	sid=$(echo $s | awk -F"_" '{print $4}')
	logger_debug "sample id: $sid"
	seqFR=($(set | grep -e "$(toupper ${NAMESPACE}_${tscs})_$sid_.*_seqfile_R" | cut -d\= -f1))	

	# sai files
	saiFR=($(ls $OUTPUT_DIR/${!s}/*.sai))
	if [[ -s  ${saiFR[0]} && -s ${saiFR[1]} ]]; then
		logger_info "The given pair of sai files does exist for the current sample ${!s}."
	else
		[[ ! -s "${saiFR[0]}" ]] && logger_warn "${saiFR[0]} file does not exist or is empty."
		[[ ! -s "${saiFR[1]}" ]] && logger_warn "${saiFR[1]} file does not exist or is empty."
		logger_warn "No pair of sai files does exist for the current sample ${!s}."
		logger_fatal "Exit the pipeline."
		exit 1;
	fi

	# error logging
	CURRENT_MAPPING_ERROR=$OUTPUT_DIR/${!s}/${!s}_mapping_err.log

	# build cli options
	bwa_sampe_cli_options=($(buildCommandLineOptions "$mapping_cmd" "$NAMESPACE" 2>$CURRENT_MAPPING_ERROR))
	rtrn=$?
	cli_options_failed_msg="[Mapping] An error occured while building the $mapping_cmd command line options for current sample ${!s}."
	exit_on_error "$CURRENT_MAPPING_ERROR" "$cli_options_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
	opts="${bwa_sampe_cli_options[@]}"
	logger_debug "[Mapping] $mapping_cmd options: $opts"
	
	# build cli
	bwa_sampe_cli="$mapping_cmd $opts ${!papamama_bwa_index_path} ${!saiFR[0]} ${!saiFR[1]} ${!seqFR[0]} ${!seqFR[1]} >$OUTPUT_DIR/${!s}/${!s}_${!ga_papamama}.sam 2>$CURRENT_MAPPING_ERROR &"

	# run the cli
	logger_debug "[Mapping] $bwa_sampe_cli"
	eval "$bwa_sampe_cli" 2>$ERROR_TMP
	pid=$!
	rtrn=$?
	eval_failed_msg="[Mapping] An error occured while eval $mapping_cmd cli."
	exit_on_error "$ERROR_TMP" "$eval_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
	logger_debug "[Mapping] $mapping_cmd pid: $pid"

	# add pid to array
	PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
done

# wait until all bwa sampe processes to finish then proceed to next step
# reinit pid array
pid_list_failed_msg="[Mapping] Failed getting process status for process $p."	
for p in "${PIDS_ARR[@]}"; do
	logger_debug "$(ps aux | grep $USER | gawk -v pid=$p '$2 ~ pid {print $0}' 2>${ERROR_TMP})"
	rtrn=$?
	exit_on_error "$ERROR_TMP" "$pid_list_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
done
logger_info "[Mapping] Wait for all $mapping_cmd processes to finish before proceed to next step."
waitalluntiltimeout "${PIDS_ARR[@]}" 2>/dev/null
logger_info "[Mapping] All $mapping_cmd processes finished. Will proceed to next step: Filtering."
PIDS_ARR=()
logger_info "[Mapping] Step completed."


#==========================================
# FILTERING
#==========================================
logger_info "[Filtering] Filter alignments."

# check sam file
logger_info "[Filtering] Check sam files output for stack samples."
for s in "${SAMPLES_STACK[@]}"; do
	logger_info "[Filtering] Current sample: ${!s}"

	samF=$(ls $OUTPUT_DIR/${!s}/*.sam)
	if [[ ! -s $samF ]]; then
		logger_fatal "[Filtering] Sam file for the current sample ${!s} does not exist or is empty."
		logger_fatal "Exit the pipeline."
		exit 1
	else
		logger_debug "[Filtering] $samF file does exist and is not empty."
	fi
done

# filter unmapped reads
logger_info "[Filtering] Filter unmapped reads."
for s in "${SAMPLES_STACK[@]}"; do
	logger_info "[Filtering] Current sample: ${!s}"

	# error logging
	CURRENT_FILTERING_ERROR=$OUTPUT_DIR/${!s}/${!s}_filtering_err.log

	samF=$(ls $OUTPUT_DIR/${!s}/*.sam)
	eval "get_mapped_reads_w_header $samF >${samF%.*}_mapped.sam 2>$CURRENT_FILTERING_ERROR &" 2>$ERROR_TMP
	pid=$!
	rtrn=$?
	eval_failed_msg="[Filtering] An error occured while eval get_mapped_reads cli."
	exit_on_error "$ERROR_TMP" "$eval_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
	logger_debug "[Filtering] get_mapped_reads pid: $pid"

	# add pid to array
	PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
done

# wait until all get_mapped_reads processes to finish then proceed to next step
# reinit pid array
pid_list_failed_msg="[Filtering] Failed getting process status for process $p."	
for p in "${PIDS_ARR[@]}"; do
	logger_debug "$(ps aux | grep $USER | gawk -v pid=$p '$2 ~ pid {print $0}' 2>${ERROR_TMP})"
	rtrn=$?
	exit_on_error "$ERROR_TMP" "$pid_list_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
done
logger_info "[Filtering] Wait for all get_mapped_reads processes to finish before proceed to next step."
waitalluntiltimeout "${PIDS_ARR[@]}" 2>/dev/null
logger_info "[Filtering] All get_mapped_reads processes finished. Will proceed to next filtering step."
PIDS_ARR=()

# check mapped sam file
logger_info "[Filtering] Check mapped sam files output for stack samples."
for s in "${SAMPLES_STACK[@]}"; do
	logger_info "[Filtering] Current sample: ${!s}"

	samF=$(ls $OUTPUT_DIR/${!s}/*_mapped.sam)
	if [[ ! -s $samF ]]; then
		logger_fatal "[Filtering] Mapped sam file for the current sample ${!s} does not exist or is empty."
		logger_fatal "Exit the pipeline."
		exit 1
	else
		logger_debug "[Filtering] $samF file does exist and is not empty."
	fi
done

# filter non unique match
logger_info "[Filtering] Filter non unique matches: extract header"
for s in "${SAMPLES_STACK[@]}"; do
	logger_info "[Filtering] Current sample: ${!s}"

	# error logging
	CURRENT_FILTERING_ERROR=$OUTPUT_DIR/${!s}/${!s}_filtering_err.log

	samF=$(ls $OUTPUT_DIR/${!s}/*_mapped.sam)
	eval "samtools view -H -S $samF >${samF%.*}_uniqMatch.sam.tmp 2>$CURRENT_FILTERING_ERROR &" 2>$ERROR_TMP
	pid=$!
	rtrn=$?
	eval_failed_msg="[Filtering] An error occured while eval get_mapped_reads cli."
	exit_on_error "$ERROR_TMP" "$eval_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
	logger_debug "[Filtering] samtools view pid: $pid"

	# add pid to array
	PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
done
logger_info "[Filtering] Wait for all samtools view processes to finish before proceed to next step."
waitalluntiltimeout "${PIDS_ARR[@]}" 2>/dev/null
logger_info "[Filtering] All samtools view processes finished. Will proceed to next filtering step."
PIDS_ARR=()

logger_info "[Filtering] Filter non unique matches: apply filter"
for s in "${SAMPLES_STACK[@]}"; do
	logger_info "[Filtering] Current sample: ${!s}"

	# error logging
	CURRENT_FILTERING_ERROR=$OUTPUT_DIR/${!s}/${!s}_filtering_err.log

	samF=$(ls $OUTPUT_DIR/${!s}/*_mapped.sam)
	eval "get_unique_matches $samF >>${samF%.*}_uniqMatch.sam.tmp 2>$CURRENT_FILTERING_ERROR &" 2>$ERROR_TMP
	pid=$!
	rtrn=$?
	eval_failed_msg="[Filtering] An error occured while eval get_unique_matches cli."
	exit_on_error "$ERROR_TMP" "$eval_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
	logger_debug "[Filtering] get_unique_matches pid: $pid"

	# add pid to array
	PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
done
logger_info "[Filtering] Wait for all get_unique_matches processes to finish before proceed to next step."
waitalluntiltimeout "${PIDS_ARR[@]}" 2>/dev/null
logger_info "[Filtering] All get_unique_matches processes finished. Will proceed to next filtering step."
PIDS_ARR=()

logger_info "[Filtering] Filter non unique matches: rename tmp to sam"
for s in "${SAMPLES_STACK[@]}"; do
	logger_info "[Filtering] Current sample: ${!s}"

	# error logging
	CURRENT_FILTERING_ERROR=$OUTPUT_DIR/${!s}/${!s}_filtering_err.log

	samF=$(ls $OUTPUT_DIR/${!s}/*_uniqMatch.sam.tmp)
	eval "mv $samF ${samF%.*} 2>$CURRENT_FILTERING_ERROR &" 2>$ERROR_TMP
	pid=$!
	rtrn=$?
	eval_failed_msg="[Filtering] An error occured while eval mv cli."
	exit_on_error "$ERROR_TMP" "$eval_failed_msg" $rtrn "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" $SESSION_TAG $EMAIL
	logger_debug "[Filtering] mv pid: $pid"

	# add pid to array
	PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
done
logger_info "[Filtering] Wait for all mv processes to finish before proceed to next step."
waitalluntiltimeout "${PIDS_ARR[@]}" 2>/dev/null
logger_info "[Filtering] All mv processes finished. Will proceed to next step: cleaning"
PIDS_ARR=()

# Clean Filtering
logger_info "[Filtering] Clean tmp files"
logger_debug "[Filtering] Remove all *.${!ga_papamama}.sam files:"
logger_debug "$(find $OUTPUT_DIR -name "*.${!ga_papamama}.sam")"
find $OUTPUT_DIR -name "*.${!ga_papamama}.sam" -delete
logger_debug "[Filtering] Remove all *_mapped.sam files:"
logger_debug "$(find $OUTPUT_DIR -name "*_mapped.sam")"
find $OUTPUT_DIR -name "*_mapped.sam" -delete
logger_info "[Filtering] All tmp files removal completed."







# close all appenders
appender_exists stderr && appender_close stderr
appender_exists console && appender_close console
appender_exists debuggerF && appender_close debuggerF























