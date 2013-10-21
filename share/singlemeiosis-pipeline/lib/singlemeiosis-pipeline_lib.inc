#! /bin/bash
#
# SINGLE MEIOSIS FUNCTIONS LIBRARY
#
# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>

# Date: 2013-10-15


#--------------------------------------------------------------
# Function: get_config_sections config_file
# Purpose : Print the list of config sections from config file
get_config_sections()
{
    typeset config_file=$1

    egrep "^\[" ${config_file} | awk -F '[][]' '
      NF==3 && $0 ~ /^\[.*\]/ { print $2 }
      ' 2>&1
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo -e "$(date '+%Y%m%d %r') [$0] Failed Error parsing configuration sections in ${config_file}." >&2
	return $rtrn
    else
	return $rtrn
    fi
}

#-------------------------------------------------------------------------
# Function : set_config_params config_file config_section [config_prefix]
# Purpose  : Set variables (prefixed by config_prefix and corresponding config_section) from config_section in config file
set_config_params()
{
    typeset config_file=$1
    typeset config_section=$2
    typeset config_prefix=$3
    typeset config_params

    [[ -n $config_prefix ]] && config_prefix=$config_prefix"_"
    config_params=$( 
       awk -F= -v Config="${config_section}" -v Prefix="${config_prefix}" '
       BEGIN { 
          Config = toupper(Config);
					Prefix = toupper(Prefix);
          patternConfig = "\\[" Config "]";
       }
       toupper($0) ~ patternConfig,(/^\[/ && toupper($0) !~ patternConfig)  { 
          if (/^#/ || /^\[/ || NF <2) next;
          sub(/^[[:space:]]*/, "");
          sub(/[[:space:]]*=[[:space:]]/, "=");
          print Prefix Config"_"$0;
       } ' ${config_file} 2>&1)
    rtrn=$?

    if [[ $rtrn -ne 0 ]]; then
	echo -e "$(date '+%Y%m%d %r') [set_config_params] Failed Error parsing configuration parameters in ${config_file}." >&2
	return $rtrn
    else
	eval "${config_params}" 2>&1
        #echo -e "${config_params[@]}"
	rtrn=$?
	if [[ $rtrn -ne 0 ]]; then
	    echo -e "$(date '+%Y%m%d %r') [set_config_params] Failed Error setting configuration parameters from ${config_file}." >&2
	    return $rtrn
	else
	    return $rtrn
	fi
    fi
}

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
