#!/usr/bin/env python

import sys
import re

import gvfParser

def main():
	summary = {}
	for line in sys.stdin:
		if re.match("^##individual-id ID=",line):
			curID=line.split("##individual-id ID=")[1].strip()
			summary[curID]={"total":{},"novel":{}}
		if not (re.match("^#",line) or re.match("^$",line)):
			gvfdata = gvfParser.gvfDataUnit(line)
			effect = gvfdata.effect()
			transkind = gvfParser.transKind(gvfdata.baseVariant())
			try:
				summary[curID]["total"][effect] += 1
				summary[curID]["total"][transkind] += 1
			except KeyError:
				summary[curID]["total"][effect] = 0
				summary[curID]["total"][transkind] = 0
			if gvfdata.isNovel():
				try:
					summary[curID]["novel"][effect] += 1
					summary[curID]["novel"][transkind] += 1
				except KeyError:
					summary[curID]["novel"][effect] = 0
					summary[curID]["novel"][transkind] = 0
	for id in summary.keys():
		formatalign = ["non_synonymous_codon",
                               "synonymous_codon",
                               "splice",
                               "other"]
		form = []
		for kind in summary[id].keys():
			for key in formatalign:
				try:
					form.append(summary[id][kind][key])
				except KeyError:
					form.append(0)
			ti = float(summary[id][kind]["transition"])
			tv = float(summary[id][kind]["transversion"])
			titv = ti/tv
			form.append(titv)
		outline =  "%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f" %tuple(form)
		print id + "\t" + outline
if __name__ == "__main__":
	main()
