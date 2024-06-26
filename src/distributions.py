''' 
	Written by Aidan Daly, DPhil candidate, University of Oxford Department of Computer Science
	
	Part of a paper tutorial on ABC parameter estimation for the Hodgkin-Huxley 
	action potential model.
'''

import math
import random
from numpy import random

''' 
	Distributions used in ABC inference for prior and kernel function specification
'''

class AbstractDistribution(object):
	def pdf(self,x):
		raise NotImplementedError
	def cdf(self,x):
		raise NotImplementedError
	def draw(self):
		raise NotImplementedError
	def getmean(self):
		return self.mean
	def getvar(self):
		return self.var
	
class Uniform(AbstractDistribution):
	def __init__(self,lo=0.0,hi=1.0):
		self.hi = hi
		self.lo = lo
		self.mean = (hi+lo)/2.0
		self.var = pow(hi-lo,2)/12.0
	def pdf(self,x):
		if x >= self.lo and x <= self.hi:
			return 1.0/(self.hi-self.lo)
		return 0.0
	def cdf(self,x):
		if x < self.lo:
			return 0.0
		if x >= self.hi:
			return 1.0
		return (x-self.lo)/(self.hi-self.lo)
	def draw(self):
		return self.lo + (self.hi-self.lo)*random.random()
	def getmean(self):
		return self.mean
	def getvar(self):
		return self.var

class Normal(AbstractDistribution):
	def __init__(self,mean=0.0,var=1.0):
		self.mean = mean
		self.var = var
	def pdf(self,x):
		return (1/(math.sqrt(2*math.pi*self.var)))*math.exp(-pow(x-self.mean,2)/(2*self.var))
	def draw(self):
		return random.normal(self.mean,math.sqrt(self.var))
	def getmean(self):
		return self.mean
	def getvar(self):
		return self.var
		
class Arbitrary(AbstractDistribution):
	def __init__(self,pool,weights=None):
		self.pool = pool
		
		if weights == None: # Defaults to uniform weights
			self.weights = [1.0/len(self.pool)]*len(self.pool)
		else: # Ensure weights are normalized
			tot_wt = sum(weights)
			self.weights = [w/tot_wt for w in weights]

		wts = self.weights
			
		mean = [0]*len(self.pool[0])
		for i,p in enumerate(self.pool):
			for j in range(len(p)):
				mean[j] = mean[j]+p[j]*wts[i]
		self.mean = mean
		
		var = [0]*len(self.pool[0])
		for i,p in enumerate(self.pool):
			for j in range(len(p)):
				var[j] = var[j]+pow(p[j]-mean[j],2.0)*wts[i]
		self.var = var
				
	def draw(self):
		r = random.rand()
		curr = 0.0
		for i,w in enumerate(self.weights):
			curr = curr+w
			if curr >= r:
				return self.pool[i]
	def getmean(self):
		return self.mean
	def getvar(self):
		return self.var
		