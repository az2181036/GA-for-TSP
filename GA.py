# -*- coding: utf-8 -*-
import random
from TSP_readfile import get_distance
from Genome import Genome

distance,length = get_distance("TSPINF.txt")
itrNum = 100  #for climbing method
INF = 0xfffffff

class GA(object): 
	"""
	种群类
	"""
	def __init__(self, popsize,chromolength,crossrate,mutationrate):
		'''
		popSize   人口总数
		chromoLength 基因长度
		crossRate 交叉概率
		mutationRare 变异概率
		population 种群列表（存放个体基因）
		totalFitness 种群总适应度
		generation 代数
		bestFiness 最佳适应度
		bestGenome 最佳个体基因
		'''
		self.popSize = popsize
		self.chromoLength = chromolength
		self.crossRate = crossrate
		self.mutationRate = mutationrate
		
		self.generation = 1
		self.population = list()
		self.totalFitness = 0 
		self.bestFitness = -1
		self.bestGenome = None

		self.InitPop()

	def InitPop(self):    #初始化种群
		for i in range(self.popSize):
			genome = list()
			for j in range(self.chromoLength-1):
				j = j + 2
				genome.append(j)
			random.shuffle(genome)
			genome.insert(0,1)
			self.population.append(Genome(genome))

	def SetPopFitness(self):	        #获取种群适应度
		self.totalFitness = 0 
		for genome in self.population:
			genome.fitness = self.SetFitness(genome.gene)
			self.totalFitness += genome.fitness
			if genome.fitness > self.bestFitness:
				self.bestFitness = genome.fitness
				self.bestGenome = genome

	def SetFitness(self,genome):   #获取个体适应度
		global distance
		fitness = 0
		for i in range(self.chromoLength):
			t = genome[i]
			k = genome[(i+1)%self.chromoLength]
			fitness += distance[t-1][t,k]
		return 1.0/fitness

	def CrossOver(self,mom,dad):
		index1 = random.randint(1,self.chromoLength-1)
		index2 = random.randint(index1,self.chromoLength-1)
		tmpBaby1 = mom[index1:index2]
		tmpBaby2 = dad[index1:index2]

		tmpDad = dad[index2:]+dad[1:index1]+tmpBaby2
		for i in range(index2-index1):
			i = i + index1
			tmpDad.remove(mom[i])
		newBaby = tmpDad[self.chromoLength-index2:]+tmpBaby1[:] + tmpDad[:self.chromoLength-index2]
		newBaby.insert(0,1)
		return newBaby
	'''

	def CrossOver(self,mom,dad):
		index = random.randint(1,self.chromoLength-3)
		tmpX = mom[index]
		newBaby = list()

		for i in range(len(dad)):
			if dad[i] == tmpX:
				tmpY = dad[i-1]
				break

		if tmpY == mom[index+1] or tmpY==mom[index-1]:
			newBaby = mom[:]
		else:
			newBaby = self.ReserveGenome(mom,index,tmpY)
		return newBaby
	'''

	def Mutate(self,baby):		
		index1 = random.randint(0,self.chromoLength-1)
		index2 = random.randint(index1,self.chromoLength-1)

		tmp = baby[index1+1:index2+1]
		tmpBaby = baby[:index1+1]+tmp[::-1]+baby[index2+1:]

		tmpFitness1 = self.SetFitness(tmpBaby)
		tmpFitness2 = self.SetFitness(baby)

		if tmpFitness1 > tmpFitness2:
			return tmpBaby
		return baby 

	def Tournament_Selection(self):
		tmp = random.sample(self.population,2)
		if tmp[0].fitness>tmp[1].fitness:
			return tmp[0]
		else:
			return tmp[1]

	def ReserveGenome(self,mom,index,tmp):
		for i in range(len(mom)):
			if mom[i]==tmp:
				if i<index:
					t  = mom[i+1:index+1]
					baby = mom[:i+1]+t+mom[index+1:]
				else:
					t = mom[index+1:i+1]
					baby = mom[:index+1]+t+mom[i+1:]
		return baby

	def ClimbMethod(self,genome,fitness):
		bestV = fitness
		bestG = genome[:]
		for i in range(itrNum):
			gen = bestG[:]
			index1 = random.randint(1,len(genome)-2)
			index2 = random.randint(1,len(genome)-2)
			gen[index1],gen[index2] = gen[index2],gen[index1]
			tmpFitness = self.SetFitness(gen)
			if tmpFitness > bestV:
				bestV = tmpFitness
				bestG = gen
		return bestG

	def ClimbMethod1(self,genome,fitness):
		bestV = fitness
		bestG = genome
		for i in range(itrNum):
			cnt = 1
			gen = bestG.gene[:]
			index1 = random.randint(0,len(genome)-3)
			tmpA = gen[index1]
			tmpB = gen[index1+1]
			valAB = distance[tmpA-1][tmpA,tmpB]
			while(True):
				index2 = random.randint(1,len(genome)-2)
				tmpC = gen[index2]
				valAC = distance[tmpA-1][tmpA,tmpC]
				
				if valAC<valAB or cnt == 20:
					break
				cnt += 1

			gen[index1+1],gen[index2] = gen[index2],gen[index1+1]
			tmpFitness = self.SetFitness(gen)
			if tmpFitness > bestV:
				bestV = tmpFitness
				bestG = Genome(gen,bestV)
		return bestG,bestV

	def ClimbMethod2(self,genome,fitness):
		global INF
		bestV = fitness
		bestG = genome
		for i in range(itrNum):
			gen = bestG.gene[:]
			index1 = random.randint(1,len(genome)-3)
			tmpA = gen[index1]
			val = INF
			dis_index = distance[tmpA-1]
			
			for i in range(self.chromoLength):
				i = i + 1
				if val >= dis_index[tmpA,i] and dis_index[tmpA,i] != 0 :
					val = dis_index[tmpA,i]
					tmpB = i
			
			for i in range(self.chromoLength):
				if gen[i] == tmpB:
					index2 = i
					break

			
			if tmpB == 1 and tmpA != 1:
				gen[index1],gen[index2+1] = gen[index2+1],gen[index1]
			elif tmpB != 1:
				gen[index1+1],gen[index2] = gen[index2],gen[index1+1]

			tmpFitness = self.SetFitness(gen)
			if tmpFitness > bestV:
				bestV = tmpFitness
				bestG = Genome(gen,bestV)
		return bestG,bestV	


	'''
	def ClimbMethodForPop(self):
		for i in range(self.popSize):
			self.population[i].gene,self.population[i].fitness = self.ClimbMethod1(self.population[i].gene,self.population[i].fitness)
			if self.population[i].fitness > self.bestFitness:
				self.bestFitness = self.population[i].fitness
				self.bestGenome = self.population[i].gene
	'''

	def Epoch(self,newPop):
		self.population = newPop               #使当前种群等于传入种群
		self.SetPopFitness()		#获取种群适应度
		#self.ClimbMethodForPop()
		newPop = list()			#生成子代容器
		self.bestGenome,self.bestFitness = self.ClimbMethod2(self.bestGenome,self.bestFitness)
						#对最优个体进行爬山，选取局部最优
		#self.bestGenome,self.bestFitness = self.ClimbMethod2(self.bestGenome,self.bestFitness)
		newPop.append(self.bestGenome)    #压入最优个体
		while(len(newPop)<self.popSize):
			mom = self.Tournament_Selection()   #锦标赛选择
			dad = self.Tournament_Selection()
			#climbMom = self.ClimbMethod1(mom.gene[:],mom.fitness)
			#climbDad = self.ClimbMethod1(dad.gene[:],dad.fitness)
			rate = random.random()

			if rate < self.crossRate:
				#baby = self.CrossOver(climbMom,climbDad)
				baby = self.CrossOver(mom.gene[:],dad.gene[:])
			else:
				#baby = climbMom
				baby = mom.gene[:]

			if rate <self.mutationRate:
				baby = self.Mutate(baby)

			newPop.append(Genome(baby))    #压入子代
		return newPop