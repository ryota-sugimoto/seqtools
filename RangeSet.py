#!/usr/bin/env python

import bisect
class RangeSet:
	def __init__(self, rangelist = []):
		self.normalizedlist = []
		for range in rangelist:
			self.insert(range[0],range[1])
					
	def __contains__(self, value):
		return bisect.bisect(self.normalizedlist, value) % 2

	def insert(self, start, end):
		if start == end:
			return
		elif start > end:
			start, end = end, start
		
		startIndex = bisect.bisect(self.normalizedlist, start)
		endIndex = bisect.bisect(self.normalizedlist, end)
		insortingValues = []
		if start not in self:
			insortingValues.append(start)
		if end not in self:
			insortingValues.append(end)
		del self.normalizedlist[startIndex: endIndex]
		for value in insortingValues:
			bisect.insort(self.normalizedlist, value)
	
	def getWidth(self):
		res = 0
		count = 0
		for l in self.normalizedlist[1::2]:
			res += l - self.normalizedlist[count]
			count += 2
		return res
	
def removedContacts(list):
	start = list[0]
	end = list[1]
	try:
		begin = list[2]
	except IndexError:
		return list
	if end == begin:
		return removedContacts([start] + list[3:])
	else:
		return list[:2] + removedContacts(list[2:])

if __name__ == "__main__":
	rangeSet = RangeSet()
	while True:
		try:
			start = int(raw_input("start: "))
			end = int(raw_input("end: "))
		except ValueError:
			print "please input integer value"
			start, end = 0, 0
		rangeSet.insert(start, end)
		print rangeSet.normalizedlist
		try:
			value = int(raw_input("value: ")) 
			if value in rangeSet:
				print value, "insides rangeSet"
			else:
				print value, "outsides rangeSet"
		except ValueError:
			pass
		print "width: ", rangeSet.getWidth()
