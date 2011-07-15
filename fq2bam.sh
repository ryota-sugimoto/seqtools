#!/usr/bin/env bash
USAGE="Usage: fq2bam.sh [options] in_seq_1 in_seq_2"

#perse options
OPTERR=0
while getopts 'a:c:o:pqur:t:l:' OPTION
do
    case $OPTION in
    o) OUT="${OPTARG%/}/" ;;
    r) REFERENCE=$OPTARG ;;
    p) CREATEPILEUP=1 ;;
    q) DO_NOT_QUALITYTRIMM=1 ;;
    u) DO_NOT_UNPAIREDFILTER=1 ;;
    a) CLIPPING_ADAPTER=$OPTARG ;;
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
checkcommand fastx_clipper

#create output directory
FILE1DIR="$(dirname "${1}")/"
[ $OUT ] || OUT="${FILE1DIR}out/"
mkdir -p "${OUT}" \
  || { echo "Failed to create output directory" >&2; exit 1; }

#path to the reference file
[ $REFERENCE ] \
  || REFERENCE="/usr/local/share/doc/hg19/hg19.fa" #default reference path
test -e "${REFERENCE}" \
  || { echo "Reference file ${REFERENCE} not found." >&2; exit 1; }

#main processes
ORIGIN1=${1}
ORIGIN2=${2}
INIT_FILE1=${1}
INIT_FILE2=${2}

#quality trimming
if [ ! $DO_NOT_QUALITYTRIMM ]
then
  [ $tOPT ] || tOPT=20 #default -t value
  [ $lOPT ] || lOPT=75 #default -l value
  NEW_FILE1="${OUT}$(basename ${INIT_FILE1}).trimmed"
  NEW_FILE2="${OUT}$(basename ${INIT_FILE2}).trimmed"
  fastq_quality_trimmer -t $tOPT -l $lOPT -i "$INIT_FILE1" \
                                          -o "$NEW_FILE1"
  fastq_quality_trimmer -t $tOPT -l $lOPT -i "$INIT_FILE2" \
                                          -o "$NEW_FILE2"
  INIT_FILE1="$NEW_FILE1"
  INIT_FILE2="$NEW_FILE2"
fi

#clipping adapter
if [ $CLIPPING_ADAPTER ]
then
  NEW_FILE1="${OUT}$(basename ${INIT_FILE1}).clipped"
  NEW_FILE2="${OUT}$(basename ${INIT_FILE2}).clipped"
  fastx_clipper -a $CLIPPING_ADAPTGER -n -v -l 75 -i "$INIT_FILE1" \
                                                  -o "$NEW_FILE1"
  fastx_clipper -a $CLIPPING_ADAPTGER -n -v -l 75 -i "$INIT_FILE2" \
                                                  -o "$NEW_FILE2"
  [[ $INIT_FILE1 != $ORIGIN1 ]] && rm $INIT_FILE1
  [[ $INIT_FILE2 != $ORIGIN2 ]] && rm $INIT_FILE2
  INIT_FILE1="$NEW_FILE1"
  INIT_FILE2="$NEW_FILE2"
fi

#filtering unpaired reads
if [ ! $DO_NOT_UNPAIREDFILTER ]
then
  NEW_FILE1="${OUT}$(basename ${INIT_FILE1}).filtered"
  NEW_FILE2="${OUT}$(basename ${INIT_FILE2}).filtered"
  fastqUnpairedFilter.py "$INIT_FILE1" \
                         "$INIT_FILE2" \
                         "$NEW_FILE1" \
                         "$NEW_FILE2"
  [[ $INIT_FILE1 != $ORIGIN1 ]] && rm $INIT_FILE1
  [[ $INIT_FILE2 != $ORIGIN2 ]] && rm $INIT_FILE2
  INIT_FILE1=$NEW_FILE1
  INIT_FILE2=$NEW_FILE2
fi

#create sai
SAIFN1="${OUT}$(basename ${INIT_FILE1}).sai"
SAIFN2="${OUT}$(basename ${INIT_FILE2}).sai" 
[ $BWATHREAD ] || BWATHREAD=4 #default threads
bwa aln -t $BWATHREAD "${REFERENCE}" \
                      "$INIT_FILE1" \
                    > "$SAIFN1"

bwa aln -t $BWATHREAD "${REFERENCE}" \
                      "$INIT_FILE2" \
                    > "$SAIFN2"

#create sam
SAMFN=$(echo $(basename $INIT_FILE1) | sed 's/[12]\([^0-9][^0-9]*\)$/\1/').sam
bwa sampe "${REFERENCE}" \
          "$SAIFN1" \
          "$SAIFN2" \
          "$INIT_FILE1" \
          "$INIT_FILE2" \
        > "${OUT}${SAMFN}"

#create bam
BAMFN="${SAMFN%.sam}.bam"
samtools import "${REFERENCE}.fai" \
                "${OUT}${SAMFN}" \
                "${OUT}${BAMFN}"

#create sorted.bam
SORTEDBAMFN="${BAMFN%.bam}.sorted"
samtools sort "${OUT}${BAMFN}" "${OUT}${SORTEDBAMFN}"

samtools index "${OUT}${SORTEDBAMFN}.bam"

#create pileup
if [ $CREATEPILEUP ]
then
  PILEUPFN=${SORTEDBAMFN%.sorted}.pileup
  samtools pileup -c -f "${REFERENCE}" "${OUT}${SORTEDBAMFN}.bam" \
                      > "${OUT}${PILEUPFN}"
  awk '( $6 >= 20 ) && ( $8 >= 20 ){ print }' "${OUT}${PILEUPFN}" \
                           > "${OUT}${PILEUPFN%.pileup}.filtered.pileup"
fi
echo "fq2bam.sh: complete"
