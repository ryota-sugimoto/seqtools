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

import time
	
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
					mytimer = time.time()
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
			return "other"
		for kind in kinds:
			if re.search(",{0,1}"+kind, ve):
				return kind
		return "other"
	
	def baseVariant(self):
		try:
			var = self["attributes"]["Variant_seq"].split(",")
			ref = self["attributes"]["Reference_seq"]
		except KeyError:
			sys.stderr.write("##WARN There is no Variant_seq\n")
			sys.stderr.write("##WARN " + str(self) + "\n")
			return None
		for seq in var:
			if ref != seq:
				return ref + seq

def transKind(trans):
	transitions = set(["AG","GA","CT","TC"])
	transversions = set(["AC","CA","GT","TG","AT","TA","GC","CG"])
	if trans.upper() in transitions:
		return "transition"
	elif trans.upper() in transversions:
		return "transversion"
	else:
		return "unknown"
				
if __name__ == "__main__":
	file = open(sys.argv[1])
	mydata = []
	for line in file:
		if not (re.match("^#",line) or re.match("^$",line)):
			mydata.append(gvfDataUnit(line))
	
	for line in mydata:
		if line.kind() == None:
			print line["attributes"]["Variant_effect"]
