import re
import RangeSet

def lineBeginWithChr(line):
	return re.match('chr*',line)

def createExonDictionary(fileobj):
	fileobj.seek(0)
	lines = fileobj.readlines()
	res = {}
	for line in lines:
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
	pileupfileobj.seek(0)
	pileuplines = pileupfileobj.readlines()
	filteredlines = []
	for line in pileuplines:
		splitedline = line.strip("\n").split("\t")
		chr = splitedline[0].strip()
		position = int(splitedline[1]) - 1 #pileup to bed format
		if positionInExon(chr, position, exonRangeSet):
			filteredlines.append(line)
	return filteredlines

if __name__ == "__main__":
	import sys
	bedFN = sys.argv[1]
	exonDict = createExonDictionary(open(bedFN,"r"))
	exonRangeSet = createExonRangeSet(exonDict)
	pileupFN = sys.argv[2]
	res = pileupExonFilter(open(pileupFN,"r"),exonRangeSet)

	for line in res:
		print line.strip("\n")
