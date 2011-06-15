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
			tmp = end
			end = start
			start = tmp

		startIndex = bisect.bisect(self.normalizedlist, start)
		endIndex = bisect.bisect(self.normalizedlist, end)
		insortingValues = []
		for l in self.normalizedlist[startIndex: endIndex]:
			insortingValues.append(l)
		if start in self:
			insortingValues.append(start)
		if end in self:
			insortingValues.append(end)
		for value in insortingValues:
			bisect.insort(self.normalizedlist, value)
		bisect.insort(self.normalizedlist, start)
		bisect.insort(self.normalizedlist, end)

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


if __name__ == '__main__':
	import sys
	sys.setrecursionlimit(1000000)
	#rangeSet = RangeSet()
	#for i in range(1000):
	#	print "start?"
	#	start = int(raw_input())
	#	print "end?"
	#	end = int(raw_input())
	#	rangeSet.insert(start, end)
	#	print rangeSet.normalizedlist
	bedlines = open("data/chr1_bed", "r").readlines()
	bedRangelist = []
	for line in bedlines:
		splitedline = line.strip("\n").split("\t")
		start = int(splitedline[0])
		end = int(splitedline[1])
		bedRangelist.append([start,end])
	print bedRangelist[:10]
	bedRangeSet = RangeSet(bedRangelist)
	print bedRangeSet.normalizedlist[:10]
	print len(bedRangeSet.normalizedlist)
	filtered = removedContacts(bedRangeSet.normalizedlist)
	print "filtered"
	print filtered[:10]
	print len(filtered)
