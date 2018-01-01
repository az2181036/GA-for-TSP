# -*- coding: utf-8 -*-
class Genome(object):
	"""docstring for Genome"""
	def __init__(self,Gene=list(),Fitness = -1):
		self.gene = Gene
		self.fitness = Fitness

	def __len__(self):
		return len(self.gene)
		