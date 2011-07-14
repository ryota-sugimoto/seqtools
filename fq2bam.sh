#!/usr/bin/env bash
USAGE="Usage: fq2bam.sh [options] in_seq_1 in_seq_2"

#perse options
OPTERR=0
while getopts 'c:o:pqur:t:l:' OPTION
do
    case $OPTION in
    o) OUT="${OPTARG%/}/" ;;
    r) REFERENCE=$OPTARG ;;
    p) CREATEPILEUP=1 ;;
    q) DO_QUALITYTRIMM=1 ;;
    u) DO_UNPAIREDFILTER=1 ;;
    t) tOPT=$OPTARG ;;
    l) lOPT=$OPTARG ;;
    c) BWATHREAD=$OPTARG ;; 
    ?) { echo $USAGE >&2 ; exit 1; };;
    esac
done
shift $(($OPTIND - 1))

#check input file existece
test -e "${1}" || { echo "File ${1} not found." >&2; exit 1; }
test -e "${2}" || { echo "File ${2} not found." >&2; exit 1; }

#check required commands
function checkcommand()
{
    which "${1}" > /dev/null 2>&1 \
      || { echo "Command ${1} not found." >&2; exit 1; }
}
checkcommand bwa
checkcommand samtools
checkcommand fastq_quality_trimmer
checkcommand fastqUnpairedFilter.py

#create output directory
FILE1DIR="$(dirname "${1}")/"
[ $OUT ] || OUT="${FILE1DIR}out/"
mkdir -p "${OUT}" \
  || { echo "Failed to create output directory" >&2; exit 1; }

FILE1=$(basename "${1}")
FILE2=$(basename "${2}")

#path to the reference file
[ $REFERENCE ] || REFERENCE="/usr/local/share/doc/hg19/hg19.fa"
test -e "${REFERENCE}" \
  || { echo "Reference file ${REFERENCE} not found." >&2; exit 1; }

#main processes
[ $tOPT ] || tOPT=20
[ $lOPT ] || lOPT=75

INIT_FILE1=${1}
INIT_FILE2=${2}

if [ $DO_QUALITYTRIMM ]
then
  fastq_quality_trimmer -t $tOPT -l $lOPT -i "$INIT_FILE1" \
                                          -o "${OUT}${FILE1}.trimmed"
  fastq_quality_trimmer -t $tOPT -l $lOPT -i "$INIT_FILE2" \
                                          -o "${OUT}${FILE2}.trimmed"
  INIT_FILE1="${OUT}${FILE1}.trimmed"
  INIT_FILE2="${OUT}${FILE2}.trimmed"
fi

if [ $DO_UNPAIREDFILTER ]
then
  fastqUnpairedFilter.py "$INIT_FILE1" \
                         "$INIT_FILE2" \
                         "${OUT}${FILE1}.filtered" \
                         "${OUT}${FILE2}.filtered"
  INIT_FILE1="${OUT}${FILE1}.filtered"
  INIT_FILE2="${OUT}${FILE2}.filtered"
fi

[ $BWATHREAD ] || BWATHREAD=4
bwa aln -t $BWATHREAD "${REFERENCE}" \
                      "$INIT_FILE1" \
                    > "${OUT}${FILE1}.sai"

bwa aln -t $BWATHREAD "${REFERENCE}" \
                      "$INIT_FILE2" \
                    > "${OUT}${FILE2}.sai"

SAMFN=$(echo $FILE1 | sed 's/[12]\([^0-9][^0-9]*\)$/\1/').sam
bwa sampe "${REFERENCE}" \
          "${OUT}${FILE1}.sai" \
          "${OUT}${FILE2}.sai" \
          "$INIT_FILE1" \
          "$INIT_FILE2" \
        > "${OUT}${SAMFN}"

BAMFN="${SAMFN%.sam}.bam"
samtools import "${REFERENCE}.fai" \
                "${OUT}${SAMFN}" \
                "${OUT}${BAMFN}"

SORTEDBAMFN="${BAMFN%.bam}.sorted"
samtools sort "${OUT}${BAMFN}" "${OUT}${SORTEDBAMFN}"

samtools index "${OUT}${SORTEDBAMFN}.bam"

if [ $CREATEPILEUP ]
then
  PILEUPFN=${SORTEDBAMFN%.sorted}.pileup
  samtools pileup -c -f "${REFERENCE}" "${OUT}${SORTEDBAMFN}.bam" \
                      > "${OUT}${PILEUPFN}"
  awk '( $6 >= 20 ) && ( $8 >= 20 ){ print }' "${OUT}${PILEUPFN}" \
                           > "${OUT}${PILEUPFN%.pileup}.filtered.pileup"
fi
echo "fq2bam.sh: complete"
