#!/usr/bin/env bash
USAGE=$(cat << EOF
Usage: fq2bam.sh [options] in_seq_1 in_seq_2

        -o dir    Specify output directory.
 
        -r file   Specify Reference fasta.

        -p        create pileup file additionally.

        -q        Do not quality trimming.

        -u        Do not unpaired filtering.

        -a        Do not adapter clipping.

        -A [ATGC][ATGC]* 
                  Spceify clipping adapter sequence.
        -M int
                  Required minimum adaptor length.
        -t int    t value of quality trimming.

        -l int    l value of quality trimming.

        -c int    Number of threads of bwa align program.

	-Q int    -Q parameter for clipper and trimmer

        -i        Illumina 1.3+ read format.

        -d        Remove internal files.
EOF
)

#default values
QUALITY_TRIM_T=20 
QUALITY_TRIM_L=70
ADAPTERSEQ="AGATCGG"
BWATHREAD=4 
OUT="$(pwd)/"
REFERENCE="/usr/local/share/doc/hg19/hg19.fa"

#perse options
while getopts 'ac:o:pqur:t:l:idA:M:Q:' OPTION
do
    case $OPTION in
    o) OUT="${OPTARG%/}/" ;;
    r) REFERENCE=$OPTARG ;;
    p) CREATEPILEUP=1 ;;
    q) DO_NOT_QUALITYTRIM=1 ;;
    u) DO_NOT_UNPAIREDFILTER=1 ;;
    a) DO_NOT_ADAPTERCLIP=1 ;;
    t) QUALITY_TRIM_T=$OPTARG ;;
    l) QUALITY_TRIM_L=$OPTARG ;;
    c) BWATHREAD=$OPTARG ;;
    i) DO_ILL2SANGER="-I" ;;
    d) RM_INTERNAL_FILES=1 ;;
    A) [ $(echo ${OPTARG} | tr -d "ATGC") ] \
         && { echo "illegal adapter sequence ${OPTARG}";  exit 1;} \
         || ADAPTERSEQ=${OPTARG} ;;
    M) MIN_ADAPTOR="-M ${OPTARG}" ;;
    Q) Q_PARAMETER="-Q ${OPTARG}" ;;
    ?) { echo -e "$USAGE" >&2 ; exit 1; };;
    esac
done
shift $(($OPTIND - 1))

#check input file existece
[ $# != 2 ] && { echo -e "$USAGE"; exit 1; }
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
 mkdir -p "${OUT}" || exit 1

#check reference file
test -e "${REFERENCE}" \
  || { echo "Reference file ${REFERENCE} not found." >&2; exit 1; }

#main processes
ORIGIN1=${1}
ORIGIN2=${2}
INIT_FILE1=${1}
INIT_FILE2=${2}

#clipping adapter
if [ ! $DO_NOT_ADAPTERCLIP ]
then
  NEW_FILE1="${OUT}$(basename ${INIT_FILE1}).clipped"
  NEW_FILE2="${OUT}$(basename ${INIT_FILE2}).clipped"
  COMMAND1=$(cat << EOF
fastx_clipper 
    -a $ADAPTERSEQ
    ${MIN_ADAPTOR}
    ${Q_PARAMETER}
    -n -v 
    -l 70 
    -i ${INIT_FILE1}
    -o ${NEW_FILE1}
EOF
)
  COMMAND2=$(cat << EOF
fastx_clipper 
    -a $ADAPTERSEQ 
    ${MIN_ADAPTOR}
    ${Q_PARAMETER}
    -n -v 
    -l 70 
    -i ${INIT_FILE2}
    -o ${NEW_FILE2}
EOF
)
  echo "${COMMAND1}"
  ${COMMAND1} || exit 1
  echo
  echo "${COMMAND2}"
  ${COMMAND2} || exit 2
  INIT_FILE1="$NEW_FILE1"
  INIT_FILE2="$NEW_FILE2"
  echo
  echo
fi

#quality trimming
if [ ! $DO_NOT_QUALITYTRIM ]
then
  NEW_FILE1="${OUT}$(basename ${INIT_FILE1}).trimmed"
  NEW_FILE2="${OUT}$(basename ${INIT_FILE2}).trimmed"
  COMMAND1=$(cat << EOF
fastq_quality_trimmer 
    -t ${QUALITY_TRIM_T}
    -v 
    ${Q_PARAMETER}
    -l ${QUALITY_TRIM_L}
    -i ${INIT_FILE1}
    -o ${NEW_FILE1}
EOF
)
  COMMAND2=$(cat << EOF
fastq_quality_trimmer 
    -t ${QUALITY_TRIM_T}
    -v 
    ${Q_PARAMETER}
    -l ${QUALITY_TRIM_L}
    -i ${INIT_FILE2}
    -o ${NEW_FILE2}
EOF
)
  echo "${COMMAND1}"
  ${COMMAND1} || exit 1
  echo
  echo "${COMMAND2}"
  ${COMMAND2} || exit 1
  [ ${RM_INTERNAL_FILES} ] && [ ${INIT_FILE1} != ${ORIGIN1} ] \
    && rm ${INIT_FILE1}
  [ ${RM_INTERNAL_FILES} ] && [ ${INIT_FILE2} != ${ORIGIN2} ] \
    && rm ${INIT_FILE2}
  INIT_FILE1="$NEW_FILE1"
  INIT_FILE2="$NEW_FILE2"
  echo
  echo
fi


#filtering unpaired reads
if [ ! $DO_NOT_UNPAIREDFILTER ]
then
  NEW_FILE1="${OUT}$(basename ${INIT_FILE1}).filtered"
  NEW_FILE2="${OUT}$(basename ${INIT_FILE2}).filtered"
  COMMAND=$(cat << EOF
fastqUnpairedFilter.py ${INIT_FILE1}
                         ${INIT_FILE2}
                         ${NEW_FILE1}
                         ${NEW_FILE2}
EOF
)
  echo "${COMMAND}"
  ${COMMAND}
  [ ${RM_INTERNAL_FILES} ] && [ ${INIT_FILE1} != ${ORIGIN1} ] \
    && rm ${INIT_FILE1}
  [ ${RM_INTERNAL_FILES} ] && [ ${INIT_FILE2} != ${ORIGIN2} ] \
    && rm ${INIT_FILE2}
  INIT_FILE1=$NEW_FILE1
  INIT_FILE2=$NEW_FILE2
  echo
  echo
fi

#convert illumina to sanger format
#if [ $DO_ILL2SANGER ]
#then
#  NEW_FILE1="${OUT}$(basename ${INIT_FILE1}).sanger"
#  NEW_FILE2="${OUT}$(basename ${INIT_FILE2}).sanger"
#  maq ill2sanger "$INIT_FILE1" \
#                 "$NEW_FILE1"
#  maq ill2sanger "$INIT_FILE2" \
#                 "$NEW_FILE2"
#  INIT_FILE1=$NEW_FILE1
#  INIT_FILE2=$NEW_FILE2
#fi


#create sai
SAIFN1="${OUT}$(basename ${INIT_FILE1}).sai"
SAIFN2="${OUT}$(basename ${INIT_FILE2}).sai" 
COMMAND1=$(cat << EOF
bwa aln ${DO_ILL2SANGER} 
  -t ${BWATHREAD}
  ${REFERENCE}
  ${INIT_FILE1}
EOF
)
COMMAND2=$(cat << EOF
bwa aln ${DO_ILL2SANGER} 
  -t ${BWATHREAD} 
  ${REFERENCE}
  ${INIT_FILE2}
EOF
)
echo "${COMMAND1}"
${COMMAND1} > "${SAIFN1}" || exit 1
echo
echo "${COMMAND2}"
${COMMAND2} > "${SAIFN2}" || exit 1
echo
echo

#create sam
SAMFN=$(echo $(basename $INIT_FILE1) \
  | sed 's/\(^s_[1-9][0-9]*_\)[12]_/\1/').sam
COMMAND=$(cat << EOF
bwa sampe ${REFERENCE}
          ${SAIFN1}
          ${SAIFN2}
          ${INIT_FILE1}
          ${INIT_FILE2}
EOF
)
echo "${COMMAND}"
${COMMAND} > "${OUT}${SAMFN}" || exit 1
echo
echo
[ ${RM_INTERNAL_FILES} ] && [ ${INIT_FILE1} != ${ORIGIN1} ] && rm ${INIT_FILE1}
[ ${RM_INTERNAL_FILES} ] && [ ${INIT_FILE2} != ${ORIGIN2} ] && rm ${INIT_FILE2}
[ ${RM_INTERNAL_FILES} ] && rm ${SAIFN1}
[ ${RM_INTERNAL_FILES} ] && rm ${SAIFN2}

#create fai
test -e "${REFERENCE}.fai" \
  || bwa index ${REFERENCE} || exit 1

#create bam
BAMFN="${SAMFN%.sam}.bam"
COMMAND=$(cat << EOF
samtools import ${REFERENCE}.fai
                ${OUT}${SAMFN}
                ${OUT}${BAMFN}
EOF
)
echo "${COMMAND}"
${COMMAND} || exit 1
echo
echo
[ ${RM_INTERNAL_FILES} ] && rm ${OUT}${SAMFN}


#create sorted.bam
SORTEDBAMFN="${BAMFN%.bam}.sorted"
COMMAND=$(cat << EOF
samtools sort 
  ${OUT}${BAMFN}
  ${OUT}${SORTEDBAMFN}
EOF
)
echo "${COMMAND}"
${COMMAND} || exit 1
echo
echo
[ ${RM_INTERNAL_FILES} ] && rm ${OUT}${BAMFN}

samtools index "${OUT}${SORTEDBAMFN}.bam" || exit 1

#create pileup
if [ $CREATEPILEUP ]
then
  PILEUPFN=${SORTEDBAMFN%.sorted}.pileup
  samtools pileup -c -f "${REFERENCE}" "${OUT}${SORTEDBAMFN}.bam" \
                      > "${OUT}${PILEUPFN}" || exit 1
  awk '( $6 >= 20 ) && ( $8 >= 20 ){ print }' "${OUT}${PILEUPFN}" \
    > "${OUT}${PILEUPFN%.pileup}.filtered.pileup" || exit 1
fi
echo "Done."
