#!/usr/bin/env python

import sys
import re

class gvfFormatError(Exception):
	def __init__(self,str):
		self.str = str
	def __str__(self):
		msg = "\n##ERROR Unexpected gvf format.\n"
		failed_str = "##ERROR input string: " + self.str + "\n"
		return msg + failed_str
		

class gvfAttributes:
	def __init__(self, in_str):
		self.data = {}
		for tag in in_str.strip("\n").split(";"):
			if tag:
				try:
					(left, right) = tag.split("=")
					self.data[left] = right
				except ValueError:
					print "error at attributes"
					raise gvfFormatError(in_str)
	
	def __getitem__(self, key):
		return self.data[key]

	def __repr__(self):
		return ";".join([ "=".join(pair) for pair in zip(self.data.keys(),self.data.values())])

class gvfDataUnit:
	def __init__(self, in_str):
		self.tags = ("chr", "source", "type", "start", "end",
                        "score", "strand", "reserved", "attributes")
		self.data = {}
		def checkData():
			start = int(self.data["start"])
			end = int(self.data["end"])
			if start > end:
				raise ValueError
			float(self.data["score"])
			availableStrand = set(["+","-",".","?"])
			if self.data["strand"] not in availableStrand:
				raise ValueError
		
		try:
			splited_line = in_str.strip("\n").split("\t")
			for (tag,value) in zip(self.tags,splited_line):
				if tag == "attributes":
					self.data[tag] = gvfAttributes(value)
				else:
					self.data[tag] = value
			checkData()
		except (ValueError,KeyError):
			raise gvfFormatError(in_str)
		
	def __getitem__(self, key):
		return self.data[key]
	
	def __repr__(self):
		return "\t".join([ str(self.data[tag]) for tag in self.tags])


	def isNovel(self):
		##This is for the human genome database.
		##You may have to edit the regular expression 
                ##for your specified gvf file.
		return not re.match("^rs",self["attributes"]["ID"])
	
	def effect(self):
		##priority is non_syn > syn > splice
		kinds = ["non_synonymous_codon", 
                         "synonymous_codon", 
                         "splice"]
		try:
			ve = self["attributes"]["Variant_effect"]
		except KeyError:
			sys.stderr.write("##WARN There is no Variant_effect\n")
			sys.stderr.write("##WARN " + str(self) + "\n")
			return
		for kind in kinds:
			if re.search(kind, ve):
				return kind
			else:
				return "other"
				
if __name__ == "__main__":
	file = open(sys.argv[1])
	mydata = []
	for line in file:
		if not (re.match("^#",line) or re.match("^$",line)):
			mydata.append(gvfDataUnit(line))
	
	for line in mydata:
		if line.kind() == None:
			print line["attributes"]["Variant_effect"]
