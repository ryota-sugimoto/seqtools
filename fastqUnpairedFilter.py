#!/usr/bin/env python

import sys

def getID(line):
	#
	#This function read a ID from the ID line in sequence data.
	#You may have to edit this function for your specified fastq format.
	#
	temp = line.split('#')[0].split(':')[-4:]
	return temp[0] + ':' + temp[1] + ':' + temp[2] + ':' + temp[3]

def getIDlist(fastqfile):
	fastqfile.seek(0)
	res = []
	count = 0
	for line in fastqfile:
		if count % 4 == 0:
			res.append(getID(line))
		count = count + 1
	return res

def fileFilter(in_fastqfile, in_fileIDlist, out_fastqfile, passSet):
	in_fastqfile.seek(0)
	out_fastqfile.seek(0)
	for ID in in_fileIDlist:
		if ID in passSet:
			for i in xrange(4):
				out_fastqfile.write(in_fastqfile.next())
		else:
			for i in xrange(4): in_fastqfile.next()

def filterMain(in_1, in_2, out_1, out_2):
	in_file1 = open(in_1,'r')
	in_file2 = open(in_2,'r')

	file1IDlist = getIDlist(in_file1)
	file2IDlist = getIDlist(in_file2)

	passSet = set(file1IDlist) & set(file2IDlist)
	
	fileFilter(in_file1,
                   file1IDlist,
                   open(out_1,"w"),
                   passSet)
	fileFilter(in_file2,
                   file2IDlist,
                   open(out_2,"w"),
                   passSet)

if __name__ == '__main__':
	usage = "Usage: fastqUnpairedFilter.py in_seq_1 in_seq_2 out_seq_1 out_seq_2"
	if not len(sys.argv) == 5:
		print usage
		exit(1)
	
	filterMain(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
