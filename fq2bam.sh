#!/bin/bash
USAGE="usage: fq2bam 1st_Seq.txt 2nd_Seq.txt [output_prefix]"
if (( $# == 0 ))
then
    echo "${USAGE}"
    exit 1
fi

if [ -e "${1}" -a -e "${2}" ]
then
    :
else
    echo "ERROR: File not exist."
    exit 1
fi

FILE1DIR="$(dirname "${1}")/"
OUT="${FILE1DIR}out/"
mkdir "${OUT}" 2> /dev/null

FILE1=$(basename "${1}")
FILE2=$(basename "${2}")

REFERENCE=~/db/hg19/hg19.fa

fastq_quality_trimmer -t 20 -l 75 -i "${1}" -o "${OUT}${FILE1}.trimmed"
fastq_quality_trimmer -t 20 -l 75 -i "${2}" -o "${OUT}${FILE2}.trimmed"

python ~/local/bin/fastqFilter.py "${OUT}${FILE1}.trimmed" \
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

if (( $# == 3 ))
then
    mv "${OUT}${BAMFN}" "${OUT}${3}.${BAMFN}"
    mv "${OUT}${SORTEDBAMFN}.bam" "${OUT}${3}.${SORTEDBAMFN}.bam"
    mv "${OUT}${PILEUPFN}" "${OUT}${3}.${PILEUPFN}"
fi

awk -f "~/local/bin/pileupFilter.awk" "${OUT}${3}.${PILEUPFN}" \
                        > "${OUT}${3}.filtered.${PILEUPFN}"

echo "fq2bam.sh: complete"
