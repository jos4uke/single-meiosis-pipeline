#! /bin/bash

#
# SINGLE MEIOSIS PIPELINE LIB FUNCTION TESTS SUITE
#

# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>

# date: 21-10-2013

#-------------------
# testFailedInArray 
#
testFailedInArray()
{
	value="vcf"
    tabix_formats=("gff" "bed" "sam" "vcf")
    if ($(in_array "$value" "${tabix_formats[@]}")); then
	echo -e "ok" >${stdoutF}
    else
	echo -e "failed" >${stderrF}
    echo -e "stderr output:"
    cat ${stderrF}
    fi

    assertTrue 'Expected output to stdout' "[ -s ${stdoutF} ]"
    assertFalse 'Unexpected output to stderr' "[ -s ${stderrF} ]"
}

#--------------------------------------
# testFailedPrintingConfigSectionsList
#
testFailedPrintingConfigSectionsList()
{
    res=$(get_config_sections $PIPELINE_DEFAULT_CONFIG 2>${stderrF})
    read -a array echo <<< $(echo -e $res | tr " " "\n")
    echo -e "config sections list: "
	for s in "${array[@]}"; do
		echo -e "- config section name: $s"
	done
    assertTrue 'Unexpected array size, should be greater than 0' "[ ${#array[@]} -gt 0 ]"
    assertFalse 'Unexpected output to stderr' "[ -s ${stderrF} ]"
}

#-----------------------------------
# testFailedReadingUserConfigParams
#
testFailedReadingUserConfigParams()
{
    for cfg in $(get_config_sections $PIPELINE_DEFAULT_CONFIG); do
	echo -e "--- Config section [${cfg}] ---"
	unset $(set | awk -F= -v cfg="${cfg}" 'BEGIN { 
          cfg = toupper(cfg);
       }
			/^cfg_/  { print $1 }') $(toupper ${cfg}_)
	set_config_params $PIPELINE_DEFAULT_CONFIG ${cfg} 2>${stderrF}
	set | grep ^$(toupper ${cfg}_)
    done	
    
    [[ -s ${stderrF} ]] && echo -e "stderr output:"; cat ${stderrF}
    assertFalse 'Unexpected output to stderr' "[ -s ${stderrF} ]"
}

#-----------------------------------
# testFailedLoadingUserConfigParams
#
testFailedLoadingUserConfigParams()
{
    prefix="varco"
    for cfg in $(get_config_sections $PIPELINE_DEFAULT_CONFIG); do
	echo -e "--- Config section [${cfg}] ---"
	unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${prefix}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)
	set_config_params $PIPELINE_DEFAULT_CONFIG ${cfg} ${prefix} 2>${stderrF}
	for params in $(set | grep ^$(toupper ${prefix}_${cfg}_)); do
	    echo -e "params: $params"
	    delimiter="="
	    declare -a arr
	    arr=($(echo ${params//$delimiter/ }))
	    echo -e "variable: ${arr[0]}"
	    echo -e "value: ${arr[1]}"
	done
	unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${prefix}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)
    done	
    
    [[ -s ${stderrF} ]] && echo -e "stderr output:"; cat ${stderrF}
    assertFalse 'Unexpected output to stderr' "[ -s ${stderrF} ]"
}


#----------------------------------
# testGetGenomeAliasesListWoSnpeff

testGetGenomeAliasesListWoSnpeff()
{
	genomes_aliases_list=($(get_genome_aliases_list_wo_snpeff $PIPELINE_DEFAULT_CONFIG 2>${stderrF} | tee ${stdoutF}))

	echo "${genomes_aliases_list[@]}" 
    for ga in "${genomes_aliases_list[@]}"
    do
	echo -e "genome alias: $ga"
    done

	assertTrue "unexpected output to standard error" "[ -n ${stderrF} ]"
	assertTrue "expected output to standard output" "[ -s ${stdoutF} ]"
	assertTrue "expected non null genome aliases list" "[ ${#genomes_aliases_list[@]} -gt 0 ]"	
}

#================

#
# Configuration
#

oneTimeSetUp()
{
    tests_start=`date +%H:%M:%S`

	. /usr/local/share/bash-common/lib/bash-common_lib.inc
    . ../share/singlemeiosis-pipeline/lib/singlemeiosis-pipeline_lib.inc
    
    TEST_OUTPUT_DIR="output"

    PIPELINE_DEFAULT_CONFIG="../share/singlemeiosis-pipeline/etc/singlemeiosis-pipeline_user.config"
    PIPELINE_USER_CONFIG="singlemeiosis-pipeline_user.config"
    TEST_SAM="data/test_mapped_MAPQ.sam"
    TEST_FASTQC_1="data/test_1_Qual_Raw_Reads_test1.fq_fastqc_summary.txt"
    TEST_FASTQC_2="data/test_2_Qual_Raw_Reads_test2.fq_fastqc_summary.txt"

	TEST_GENOMES_PATH="/data/SEQUENCES/GENOME"
	TEST_BWA_INDEXES="/data/SEQUENCES/INDEX/bwa"
	TEST_SAMTOOLS_INDEXES="/data/SEQUENCES/INDEX/samtools"
	TEST_SNPEFF_VERSION="/usr/local/src/snpEff_3_0"
	TEST_SNPEFF_DATA=$TEST_SNPEFF_VERSION/"data"
}

setUp()
{
    OUTPUT_DIR="${SHUNIT_TMPDIR}/OUTPUT"
    stdoutF="${OUTPUT_DIR}/stdoutF"
    stderrF="${OUTPUT_DIR}/stderrF"
    mkdir $OUTPUT_DIR
}

tearDown()
{  
    echo "Test starts ${tests_start}"
    tests_end=`date  +%H:%M:%S`
    echo "Test ends ${tests_end}"
    exec_start_time=`date +%s -d ${tests_start}`
    exec_end_time=`date +%s -d ${tests_end}`
    exec_time=$[${exec_end_time}-${exec_start_time}]
    echo |awk -v time="$exec_time" '{print "execution time: " strftime("%Hh:%Mm:%Ss", time, 1)}'
    rm -rf $OUTPUT_DIR
    echo "------------"
}

. /usr/local/share/shunit2-2.1.6/src/shunit2
