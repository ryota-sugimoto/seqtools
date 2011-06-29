#!/usr/bin/env bash
USAGE="Usage: pileupExonFilter.sh bedfile in_pileup out_pileup > report"
if (( $# == 0 ))
then
    echo "${USAGE}" >&2
    exit 1
fi

#check input file existece
test -e "${1}" || { echo "${1} not found." >&2; exit 1; }
test -e "${2}" || { echo "${2} not found." >&2; exit 1; }

#check required commands
function checkcommand()
{
    which "${1}" > /dev/null 2>&1 \
      || { echo "Command ${1} not found." >&2; exit 1; }
}
checkcommand pileupExonFilter.py

#main process
pileupExonFilter.py "${1}" "${2}" > "${3}" \
 || { echo "pileupExonFilter.py failed" >&2; exit 1; }

#report
for chr in `awk '{ print $1 }' < "${1}.width" | sort -n -t r -k 2 | uniq`
do
    BEDWIDTH=$(grep "${chr}[[:space:]]" "${1}.width" | awk '{print $2}')
    if [ $chr != "tortal" ]
    then
        PILEUPWIDTH=$(grep "${chr}[[:space:]]" "${3}" | wc -l)
    else
        PILEUPWIDTH=$(wc -l < "${3}")
    fi
    printf "%25s%15s%15s%15s\n" \
    $chr \
    $BEDWIDTH \
    $PILEUPWIDTH \
    $(echo "scale=6; $PILEUPWIDTH / $BEDWIDTH" | bc)
done
