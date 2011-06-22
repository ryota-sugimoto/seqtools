#!/usr/bin/env python

import sys

def getID(s):
	temp = s.split('#')[0].split(':')[-4:]
	return temp[0] + ':' + temp[1] + ':' + temp[2] + ':' + temp[3]

def getIDlist(filelines):
	res = []
	count = 0
	for line in filelines:
		if count % 4 == 0:
			res.append(getID(line))
		count = count + 1
	return res

def fileFilter(srcfilelines, fileIDlist, passSet):
	srciter = iter(srcfilelines)
	res = []
	for ID in fileIDlist:
		if ID in passSet:
			for i in xrange(4): res.append(srciter.next())
		else:
			for i in xrange(4): srciter.next()
	return res

def filterMain(src1, src2, tar1, tar2):
	srcfile1list = open(src1,'r').readlines()
	srcfile2list = open(src2,'r').readlines()

	file1IDlist = getIDlist(srcfile1list)
	file2IDlist = getIDlist(srcfile2list)

	passSet = set(file1IDlist) & set(file2IDlist)
	
	tarfile1list = fileFilter(srcfile1list, file1IDlist, passSet)
	tarfile2list = fileFilter(srcfile2list, file2IDlist, passSet)

	open(tar1,'w').writelines(tarfile1list)
	open(tar2,'w').writelines(tarfile2list)

if __name__ == '__main__':
	usage = "Usage: fastqUnpairedFilter.py in_seq_1 in_seq_2 out_seq_1 out_seq_2"
	if not len(sys.argv) == 5:
		print usage
		exit(1)
	
	filterMain(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
