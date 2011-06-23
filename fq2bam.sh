#!/usr/bin/env bash
USAGE="usage: fq2bam.sh in_seq_1 in_seq_2"
if (( $# == 0 ))
then
    echo "${USAGE}"
    exit 1
fi

#check input file existece
test -e "${1}" || { echo "${1} not found."; exit 1; }
test -e "${2}" || { echo "${2} not found."; exit 1; }

#check required commands
function checkcommand()
{
  which "${1}" > /dev/null 2>&1 || { echo "Command ${1} not found." ; exit 1; }
}
checkcommand bwa
checkcommand samtools
checkcommand fastq_quality_trimmer
checkcommand fastqUnpairedFilter.py
checkcommand pileupFilter.awk

#create output directory
FILE1DIR="$(dirname "${1}")/"
OUT="${FILE1DIR}out/"
mkdir "${OUT}" 2> /dev/null

FILE1=$(basename "${1}")
FILE2=$(basename "${2}")

#path to the reference file
REFERENCE=/usr/local/share/doc/hg19/hg19.fa
test -e "${REFERENCE}" || { echo "${REFERENCE} not found."; exit 1; }

#main processes
fastq_quality_trimmer -t 20 -l 75 -i "${1}" -o "${OUT}${FILE1}.trimmed"
fastq_quality_trimmer -t 20 -l 75 -i "${2}" -o "${OUT}${FILE2}.trimmed"

fastqUnpairedFilter.py "${OUT}${FILE1}.trimmed" \
                                  "${OUT}${FILE2}.trimmed" \
                                  "${OUT}${FILE1}.trimmed.filtered" \
                                  "${OUT}${FILE2}.trimmed.filtered"

bwa aln -t 4 "${REFERENCE}" \
             "${OUT}${FILE1}.trimmed.filtered" \
           > "${OUT}${FILE1}.sai"

bwa aln -t 4 "${REFERENCE}" \
             "${OUT}${FILE2}.trimmed.filtered" \
           > "${OUT}${FILE2}.sai"

SAMFN="${FILE1/_?_sequence/_sequence}.sam"
bwa sampe "${REFERENCE}" \
          "${OUT}${FILE1}.sai" \
          "${OUT}${FILE2}.sai" \
          "${OUT}${FILE1}.trimmed.filtered" \
          "${OUT}${FILE2}.trimmed.filtered" \
        > "${OUT}${SAMFN}"

BAMFN="${FILE1/_?_sequence.txt/_sequence}.bam"
samtools import "${REFERENCE}.fai" \
                "${OUT}${SAMFN}" \
                "${OUT}${BAMFN}"

SORTEDBAMFN="${BAMFN/.bam/.sorted}"
samtools sort "${OUT}${BAMFN}" "${OUT}${SORTEDBAMFN}"

samtools index "${OUT}${SORTEDBAMFN}.bam"

PILEUPFN="${FILE1/_?_sequence.txt/_sequence}.pileup"
samtools pileup -cf "${REFERENCE}" "${OUT}${SORTEDBAMFN}.bam" > "${OUT}${PILEUPFN}"

awk '( $6 >= 20 ) && ( $8 >= 20 ){ print }' "${OUT}${PILEUPFN}" \
                                          > "${OUT}filtered.${PILEUPFN}"

echo "fq2bam.sh: complete"
