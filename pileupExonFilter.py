#!/usr/bin/env python

import re
import RangeSet

def lineBeginWithChr(line):
	return re.match('chr*',line)

def createExonDictionary(fileobj):
	fileobj.seek(0)
	res = {}
	for line in fileobj:
		if lineBeginWithChr(line):
			splitedline = line.strip("\n").split("\t")
			try:
				chr = splitedline[0].strip()
				exonbegin = int(splitedline[1])
				exonend = int(splitedline[2]) + 1
			except (IndexError, ValueError):
				print "ERROR: Unexpected format in bed."
				print line
				raise
			try:
				res[chr].append([exonbegin, exonend])
			except KeyError:
				res[chr] = [[exonbegin, exonend]]
	return res

def createExonRangeSet(exonDict):
	exonRangeSet = {}
	for chr in exonDict.keys():
		exonRangeSet[chr] = RangeSet.RangeSet(exonDict[chr])
	return exonRangeSet
	

def positionInExon(chr, position , exonRangeSet):
	try:
		return position in exonRangeSet[chr]
	except KeyError:
		return False

def pileupExonFilter(pileupfileobj, exonRangeSet):
	import sys
	pileupfileobj.seek(0)
	for line in pileupfileobj:
		splitedline = line.strip("\n").split("\t")
		try:
			chr = splitedline[0].strip()
			position = int(splitedline[1]) - 1 #pileup2bed position
		except (IndexError, ValueError):
			print "ERROR: Unexpected format in pileup."
			print line
			raise
		if positionInExon(chr, position, exonRangeSet):
			sys.stdout.write(line)

if __name__ == "__main__":
	import sys
	usage = "Usage: pileupExonFilter.py bedfile in_pileup > out_pileup"
	if len(sys.argv) != 3:
		print usage
		exit(1)
	bedFN = sys.argv[1]
	exonDict = createExonDictionary(open(bedFN,"r"))
	exonRangeSet = createExonRangeSet(exonDict)
	pileupFN = sys.argv[2]
	pileupExonFilter(open(pileupFN,"r"),exonRangeSet)
