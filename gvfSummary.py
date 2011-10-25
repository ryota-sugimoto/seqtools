#!/usr/bin/env python

import sys
import re

import gvfParser

def main():
	summary = {}
	for line in sys.stdin:
		if re.match("^##individual-id ID=",line):
			currentID=line.split("##individual-id ID=")[0]
			
			summary[currentID]={"total":{"non_synonymous_codon":0,
                                                     "synonymous_codon":0,
                                                     "splice":0,
                                                     "other":0},
                                            "novel":{"non_synonymous_codon":0,
                                                     "synonymous_codon":0,
                                                     "splice":0,
                                                     "other":0}}
		if not (re.match("^#",line) or re.match("^$",line)):
			gvfdata = gvfParser.gvfDataUnit(line)
			summary[currentID]["total"][gvfdata.effect()] += 1
			if gvfdata.isNovel():
				summary[currentID]["novel"][gvfdata.effect()] += 1

if __name__ == "__main__":
	main()
