#-*-coding:UTF-8-*-
import os
import math

class Point(object):
	"""docstring for Point"""
	def __init__(self, t):	
		self.flag = int(t[0])
		self.x = int(t[1])
		self.y = int(t[2])

def get_inf(filename):	
	fin = open(filename)
	for line in fin:
		word = line.split()
		if 'NODE_COORD_SECTION' in word:
			break
	inf = list()
	for line in fin:
		pst = line.split()
		if len(pst)==3:
			inf.append(Point(pst))
	return inf
	
def get_distance(filename):
	inf = get_inf(filename)
	length = len(inf)
	distance = list()

	for i in range(length):
		distance_i=dict()
		for j in range(length):
			'''
			get distance_i[i][j]  j = 1...length
			'''
			A = inf[i].x,inf[i].y
			B = inf[j].x,inf[j].y
			distance_i[i+1,j+1] = distance_AB(A,B)
		distance.append(distance_i)
	return distance,length

def distance_AB(A,B):
	x = float(A[0])-float(B[0])
	y = float(A[1])-float(B[1])
	return math.sqrt(x**2+y**2)

def GetPathValue(path):
	val = 0
	for i in range(len(path)-1):
		val += distance[path[i]-1][path[i],path[i+1]]
	return val

if __name__ == '__main__':
	distance,length=get_distance("TSPINF.txt")