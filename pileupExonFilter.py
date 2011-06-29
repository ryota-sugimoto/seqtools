#!/usr/bin/env python

import re
import sys
import os.path
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
				exonbegin = int(splitedline[1]) + 1
				exonend = int(splitedline[2]) + 1
			except (IndexError, ValueError):
				sys.stderr.write("ERROR: Unexpected format in bed.\n")
				sys.stderr.write(line)
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
	pileupfileobj.seek(0)
	for line in pileupfileobj:
		splitedline = line.strip("\n").split("\t")
		try:
			chr = splitedline[0].strip()
			position = int(splitedline[1])
		except (IndexError, ValueError):
			sys.stderr.write("ERROR: Unexpected format in pileup.\n")
			sys.stderr.write(line)
			raise
		if positionInExon(chr, position, exonRangeSet):
			sys.stdout.write(line)

def writeExonWidth(bedFN, exonRangeSet):
	dir = os.path.dirname(bedFN)
	raw_bedFN = os.path.basename(bedFN)
	if dir:
		dir = dir + "/"
	else:
		dir = "./"
	widthfile = open(dir + raw_bedFN + ".width", "w")
	sum = 0
	for chr in exonRangeSet.keys():
		width = int(exonRangeSet[chr].getWidth())
		sum += width
		print >> widthfile, chr, "\t", width
	print >> widthfile, "tortal", sum
	widthfile.close()
	
def readBedfile(bedFN):
	exonDict = createExonDictionary(open(bedFN, "r"))
	exonRangeSet = createExonRangeSet(exonDict)
	writeExonWidth(bedFN, exonRangeSet)
	return exonRangeSet
		
if __name__ == "__main__":
	usage = "Usage: pileupExonFilter.py bedfile in_pileup > out_pileup"
	if len(sys.argv) != 3:
		print usage
		exit(1)
	bedFN = sys.argv[1]
	exonRangeSet = readBedfile(bedFN)
	pileupFN = sys.argv[2]
	pileupExonFilter(open(pileupFN,"r"),exonRangeSet)
