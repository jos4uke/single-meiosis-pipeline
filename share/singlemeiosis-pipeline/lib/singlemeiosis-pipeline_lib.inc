#! /bin/bash
#
# SINGLE MEIOSIS FUNCTIONS LIBRARY
#
# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>

# Date: 2013-10-15

# VERSION="RC1"

###
### Genomes Listing
###

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

###
### VCF utils
###

#----------------------------------------------
# reformat_vcf function
#
reformat_vcf()
{
	local vcf_input=$1
	awk '{print $1 ";" $2 ";" $4 ";" $8}' $vcf_input | awk -F ";|=|," '{print $1 "\t" $2 "\t" $3 "\t" $7+$8 "\t" $11 "\t" $15}'
}
