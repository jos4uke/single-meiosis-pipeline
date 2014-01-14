#! /bin/bash
#
# SINGLE MEIOSIS FUNCTIONS LIBRARY
#
# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>

# Date: 2013-10-15


#---------------------------------------------------------
# Function: get_genome_aliases_list_wo_snpeff config_file
# Purpose: Get the genome aliases list having one fasta file and bwa/samtools indexes 
#
get_genome_aliases_list_wo_snpeff()
{
	local PIPELINE_USER_CONFIG_FILE=$1
	
	local GENOMES_BASE_PATH=$(grep -e "^GENOMES_BASE_PATH=" $PIPELINE_USER_CONFIG_FILE | awk -F"=" '{print $2}')
	local INDEXES_BASE_PATH=$(grep -e "^INDEXES_BASE_PATH=" $PIPELINE_USER_CONFIG_FILE | awk -F"=" '{print $2}')
	eval "local BWA_INDEXES=$(grep -e "^BWA_INDEXES=" $PIPELINE_USER_CONFIG_FILE | awk -F"=" '{print $2}')"
	eval "local SAMTOOLS_INDEXES=$(grep -e "^SAMTOOLS_INDEXES=" $PIPELINE_USER_CONFIG_FILE | awk -F"=" '{print $2}')"

	local genomes_list_fasta=($(get_genomes_dir_list_with_one_fasta ${GENOMES_BASE_PATH} 2>/dev/null))
	local genomes_w_bwa_idx=($(get_genomes_list_with_bwa_index ${BWA_INDEXES}/$(get_tool_version "bwa") "${genomes_list_fasta[@]}" 2>/dev/null))
	local GENOME_ALIASES_LIST=($(get_genomes_list_with_samtools_index ${SAMTOOLS_INDEXES}/$(get_tool_version "samtools") "${genomes_w_bwa_idx[@]}" 2>/dev/null))

	echo "${GENOME_ALIASES_LIST[@]}"
}

#---------------------------------------------------------
# Function: get_genome_aliases_list_w_bwa config_file
# Purpose: Get the genome aliases list having one fasta file and bwa index 
#
get_genome_aliases_list_w_bwa()
{
	local PIPELINE_USER_CONFIG_FILE=$1
	
	local GENOMES_BASE_PATH=$(grep -e "^GENOMES_BASE_PATH=" $PIPELINE_USER_CONFIG_FILE | awk -F"=" '{print $2}')
	local INDEXES_BASE_PATH=$(grep -e "^INDEXES_BASE_PATH=" $PIPELINE_USER_CONFIG_FILE | awk -F"=" '{print $2}')
	eval "local BWA_INDEXES=$(grep -e "^BWA_INDEXES=" $PIPELINE_USER_CONFIG_FILE | awk -F"=" '{print $2}')"

	local genomes_list_fasta=($(get_genomes_dir_list_with_one_fasta ${GENOMES_BASE_PATH} 2>/dev/null))
	local GENOME_ALIASES_LIST=($(get_genomes_list_with_bwa_index ${BWA_INDEXES}/$(get_tool_version "bwa") "${genomes_list_fasta[@]}" 2>/dev/null))

	echo "${GENOME_ALIASES_LIST[@]}"
}

#-------------------------------------------------------------------------
# Function : exit_on_error err msg status log 
# Purpose  : check if current process has to exit when checking for exit status, if different of 0 exits with error msg send to log with err file reference for details
exit_on_error()
{
	err=$1
	msg=$2
	status=$3
	log=$4
	tag=$5
	recipient=$6
	if [[ "$status" -ne 0 ]]; then
		echo -e "$(date '+%Y_%m_%d %T') $msg" | tee -a $err 2>&1 | [[ -n $log ]] && tee -a $log 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code $status." | tee -a $err 2>&1 | [[ -n $log ]] && tee -a $log 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Pipeline error] Tail error output:" | [[ -n $log ]] && tee -a $log 2>&1
		echo -e "... $(tail -n 10 $err)" 2>&1 | [[ -n $log ]] && tee -a $log 2>&1
    	echo -e "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $err." | [[ -n $log ]] && tee -a $log 2>&1
		# send an email
		proc_error_msg="$tag process error: $msg\nExit status: $status\nMore details can be found in $(readlink -f $err)\nLast error output:\n... $(tail -n 10 $err)"
		[[ -n $recipient ]] && sendEmail $recipient "[$(basename ${0%.*})] $tag process error" "$proc_error_msg" 
		exit $status
	fi
}
