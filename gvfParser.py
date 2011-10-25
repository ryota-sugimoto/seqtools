#!/usr/bin/env python

import sys
import re

class gvfVariantEffect:
	pass

class gvfAttributes:
	def __init__(self, in_str):
		self.data = {}
		for tag in in_str.strip("\n").split(";"):
			if tag:
				try:
					(left, right) = tag.split("=")
					self.data[left] = right
				except ValueError:
					print >> sys.stderr,("##WARN",
                                                   "Unexpected ",
                                                   "format in gvf attribute: ",
                                                   in_str)
	
	def __getitem__(self, key):
		return self.data[key]

	def __repr__(self):
		return ";".join([ "=".join(pair) for pair in zip(self.data.keys(),self.data.values())])

class gvfDataUnit:
	def __init__(self, in_str):
		self.tags = ("chr", "source", "type", "start", "end",
                        "score", "strand", "reserved", "attributes")
		self.data = {}
		for (tag,value) in zip(self.tags,in_str.strip("\n").split("\t")):
			if tag == "attributes":
				self.data[tag] = gvfAttributes(value)
			else:
				self.data[tag] = value

	def isNovel(self):
		##This is for the human genome database.
		##You may have to edit the regular expression 
                ##for your specified gvf file.
		return not re.match("^rs",self["attributes"]["ID"])
	
	def kind(self):
		##priority is non_syn > syn > splice
		kinds = ["non_synonymou_codon", 
                         "synonymous_codon", 
                         "splice"]
		try:
			ve = self["attributes"]["Variant_effect"]
		except KeyError:
			print >> sys.stderr,"##WARN there is no Variant_effect"
			print >> sys.stderr,"##WARN",self
			return
		for kind in kinds:
			if re.search(kind, ve):
				return kind
				
	def __getitem__(self, key):
		return self.data[key]
	
	def __repr__(self):
		return "\t".join([ str(self.data[tag]) for tag in self.tags])

if __name__ == "__main__":
	file = open(sys.argv[1])
	mydata = []
	for line in file:
		if not (re.match("^#",line) or re.match("^$",line)):
			mydata.append(gvfDataUnit(line))
	
	for line in mydata:
		print line.kind()
