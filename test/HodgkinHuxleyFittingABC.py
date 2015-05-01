''' 
	Written by Aidan Daly, DPhil candidate, University of Oxford Department of Computer Science
	
	Part of a paper tutorial on ABC parameter estimation for the Hodgkin-Huxley 
	action potential model.
'''

import fitting
import distributions as Dist
import math

# Functional curation imports
try:
    import unittest2 as unittest
except ImportError:
    import unittest

import fc
import fc.utility.test_support as TestSupport
import fc.simulations.model as Model
from fc.simulations.solvers import CvodeSolver
from fc.language import values as Val

'''
	Fitting of the Potassium and Sodium conductance components of the full Hodgkin-Huxley
	model to time course data from all voltage clamp experiments in tandem. 
	(Voltage independent)
'''

class ABCFitting(unittest.TestCase):

	# Fits the 5 parameters controlling potassium conductance in the full Hodgkin-Huxley model
	def TestGK_fitting(self):
		proto = fc.Protocol('projects/HodgkinHuxleyABC/test/protocols/VoltageClamp.txt')
		proto.SetOutputFolder('HodgkinHuxleyABC')
		proto.SetModel('projects/HodgkinHuxleyABC/hodgkin_huxley.cellml', useNumba=False)
		proto.model.SetSolver(CvodeSolver())

		outfile = open('ABCPredPotassium.txt','w')

		priors = [Dist.Uniform(0,1),Dist.Uniform(0,100),Dist.Uniform(1,100),Dist.Uniform(0,1),Dist.Uniform(1,100)]
		init = [0.5,50.0,50.0,0.5,50.0]

		def kern(orig,new=None):
			g1 = Dist.Normal(0.0,0.1)
			g10 = Dist.Normal(0.0,1.0)
			g100 = Dist.Normal(0.0,10.0)
			
			if new==None:
				perturb = [g1.draw(),g100.draw(),g100.draw(),g1.draw(),g100.draw()]
				new = []
				for i in range(len(orig)):
					new = new + [orig[i]+perturb[i]]
				return new
			else:
				prob = 1.0
				prob = prob*g1.pdf(new[0]-orig[0])
				prob = prob*g100.pdf(new[1]-orig[1])
				prob = prob*g100.pdf(new[2]-orig[2])
				prob = prob*g1.pdf(new[3]-orig[3])
				prob = prob*g100.pdf(new[4]-orig[4])
				return prob

		def distance(params,times,vals):
			# Reset model on each call
			proto.model.ResetState()
			proto.SetInput('ank1',Val.Simple(params[0]))
			proto.SetInput('ank2',Val.Simple(params[1]))
			proto.SetInput('ank3',Val.Simple(params[2]))
			proto.SetInput('bnk1',Val.Simple(params[3]))
			proto.SetInput('bnk2',Val.Simple(params[4]))
			proto.Run(verbose=False,writeOut=False)
			proto_time = proto.outputEnv.LookUp('time').array
			proto_G_K = proto.outputEnv.LookUp('G_K').array

			# Iterate through voltage clamp experiments and average RMSE
			tot_rmse = 0.0
			for i in range(11):
				tot_rmse = tot_rmse + CheckAgainst(proto_time[i],proto_G_K[i],times[i],vals[i])
			return tot_rmse/11

		# Compose a 2D array of the original Hodgkin-Huxley time and conductance data
		original_times = []
		original_G_k = []
		for i in range(11):
			[times,vals,a,b] = HodgkinHuxley.fig3(i)
			original_times = original_times + [times]
			original_G_k = original_G_k + [vals]

		result = fitting.approx_bayes_smc_adaptive(init,priors,original_times,original_G_k,prior_func,kern,distance,100,10000,0.003)
		
		# Write results to standard output and ABCPredPotassium.txt
		print result.getmean()
		print result.getvar()
		outfile.write(str(result.pool)+"\n")
		outfile.write(str(result.getmean())+"\n")
		outfile.write(str(result.getvar())+"\n")

	# Fits the 9 parameters controlling Sodium conductance in the full Hodgkin-Huxley model
	def TestGNa_fitting(self):
		proto = fc.Protocol('projects/HodgkinHuxleyABC/test/protocols/VoltageClamp.txt')
		proto.SetOutputFolder('HodgkinHuxleyABC')
		proto.SetModel('projects/HodgkinHuxleyABC/hodgkin_huxley.cellml', useNumba=False)
		proto.model.SetSolver(CvodeSolver())

		outfile = open('ABCPredSodium.txt','w')

		priors_m = [Dist.Uniform(0,1),Dist.Uniform(0,100),Dist.Uniform(1,100),Dist.Uniform(0,10),Dist.Uniform(1,100)]
		priors_h = [Dist.Uniform(0,1),Dist.Uniform(1,100),Dist.Uniform(0,100),Dist.Uniform(1,100)]
		priors = priors_m + priors_h
		init = [0.5,50.0,50.0,0.5,50.0,0.5,50.0,50.0,50.0]

		def kern(orig,new=None):
			g1 = Dist.Normal(0.0,0.1)
			g10 = Dist.Normal(0.0,1.0)
			g100 = Dist.Normal(0.0,10.0)
			
			if new==None:
				perturb_m = [g1.draw(),g100.draw(),g100.draw(),g1.draw(),g100.draw()]
				perturb_h = [g1.draw(),g100.draw(),g100.draw(),g100.draw()]
				perturb = perturb_m + perturb_h
				new = []
				for i in range(len(orig)):
					new = new + [orig[i]+perturb[i]]
				return new
			else:
				prob = 1.0
				prob = prob*g1.pdf(new[0]-orig[0])
				prob = prob*g100.pdf(new[1]-orig[1])
				prob = prob*g100.pdf(new[2]-orig[2])
				prob = prob*g1.pdf(new[3]-orig[3])
				prob = prob*g100.pdf(new[4]-orig[4])
				prob = prob*g1.pdf(new[5]-orig[5])
				prob = prob*g100.pdf(new[6]-orig[6])
				prob = prob*g100.pdf(new[7]-orig[7])
				prob = prob*g100.pdf(new[8]-orig[8])
				return prob

		def distance(params,times,vals):
			# Reset model on each call
			proto.model.ResetState()
			proto.SetInput('amk1',Val.Simple(params[0]))
			proto.SetInput('amk2',Val.Simple(params[1]))
			proto.SetInput('amk3',Val.Simple(params[2]))
			proto.SetInput('bmk1',Val.Simple(params[3]))
			proto.SetInput('bmk2',Val.Simple(params[4]))
			proto.SetInput('ahk1',Val.Simple(params[5]))
			proto.SetInput('ahk2',Val.Simple(params[6]))
			proto.SetInput('bhk1',Val.Simple(params[7]))
			proto.SetInput('bhk2',Val.Simple(params[8]))
			proto.Run(verbose=False,writeOut=False)
			proto_time = proto.outputEnv.LookUp('time').array
			proto_G_Na = proto.outputEnv.LookUp('G_Na').array

			# Iterate through voltage clamp experiments and average RMSE
			tot_rmse = 0.0
			for i in range(12):
				tot_rmse = tot_rmse + CheckAgainst(proto_time[i],proto_G_Na[i],times[i],vals[i])
			return tot_rmse/12

		# Compose a 2D array of the original Hodgkin-Huxley time and conductance data
		original_times = []
		original_G_na = []
		for i in range(12):
			[times,vals,a,b,c,d] = HodgkinHuxley.fig6(i)
			original_times = original_times + [times]
			original_G_na = original_G_na + [vals]

		result = fitting.approx_bayes_smc_adaptive(init,priors,original_times,original_G_na,prior_func,kern,distance,100,10000,0.003)
		
		# Write results to standard output and ABCPredSodium.txt
		print result.getmean()
		print result.getvar()
		outfile.write(str(result.pool)+"\n")
		outfile.write(str(result.getmean())+"\n")
		outfile.write(str(result.getvar())+"\n")

'''
	HELPER METHODS
'''

# Evaluates RMSE between experimental and predicted values
# Uses time points in simulation that are closest to experimental
def CheckAgainst(predTimes, predVals, experTimes, experVals):
    curr = 0
    predValsClosest = []
    predTimesClosest = []
        
    # Finds experimental output at times closest to experimental
    for tval in experTimes:
        while tval > predTimes[curr+1]:
            if curr >= len(predTimes)-2:
                break
            curr = curr+1
        if abs(tval-predTimes[curr]) < abs(tval-predTimes[curr+1]):
        	#print predTimes[curr]
        	predTimesClosest = predTimesClosest + [predTimes[curr]]
        	predValsClosest = predValsClosest + [predVals[curr]]
        else:
        	#print predTimes[curr+1]
        	predTimesClosest = predTimesClosest + [predTimes[curr+1]]
        	predValsClosest = predValsClosest + [predVals[curr+1]]
    # Calculate squared error
    sq_err = 0

    for i,val in enumerate(experVals):
        sq_err = sq_err + pow(val-predValsClosest[i],2)
    return math.sqrt(sq_err/len(experVals))

# Simple multiplicative prior for a list of independent Distribution objects
def prior_func(priors,params):
	prob = 1.0
	for i,distr in enumerate(priors):
		prob = prob * distr.pdf(params[i])		
	return prob