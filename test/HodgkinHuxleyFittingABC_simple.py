''' 
    Written by Aidan Daly, DPhil candidate, University of Oxford Department of Computer Science
    
    Part of a paper tutorial on ABC parameter estimation for the Hodgkin-Huxley 
    action potential model.
'''

import fitting
import distributions

try:
    import unittest2 as unittest
except ImportError:
    import unittest

import fc
import fc.utility.test_support as TestSupport
import fc.simulations.model as Model
from fc.simulations.solvers import CvodeSolver
from fc.language import values as V

'''
    Fitting of the Potassium and Sodium conductance components of the simplified Hodgkin-Huxley
    model to time course data from each of a series of voltage clamp experiments.
    (Voltage dependent) 
'''

class TestHHProto(unittest.TestCase):
            
    # Fits Potassium conductance gating rates alpha_n and beta_n from the simplified 
    #   Hodgkin-Huxley model over a range of depolarization values.
    def TestGKFitting(self):
        proto = fc.Protocol('projects/HodgkinHuxleyABC/test/protocols/VoltageClamp_simple.txt')
        proto.SetOutputFolder('HodgkinHuxleyABC_alphabeta')
        proto.SetModel('projects/HodgkinHuxleyABC/hodgkin_huxley.cellml', useNumba=False)
        proto.model.SetSolver(CvodeSolver())

        init = [0.5,0.5]
        priors = [distributions.Uniform()]*2

        # ABC expects this form - sets alpha/beta, runs protocol, then returns sq_err of result
        def distance(params,times,vals):
            # Need to reset model state before evaluating again
            proto.model.ResetState()
            proto.SetInput('alpha_n',V.Simple(params[0]))
            proto.SetInput('beta_n',V.Simple(params[1]))
            proto.Run(verbose=False,writeOut=False)
            return CheckAgainst(proto.outputEnv.LookUp('time').array,proto.outputEnv.LookUp('G_K').array,times,vals)
        
        # Kernel function is random walk with spread of 0.1
        def kern(orig,new=None):
            gauss = distributions.Normal(0.0,0.1)
            if new == None:
                new = [param + gauss.draw() for param in orig]
                return new
            else:
                prob = 1.0
                for i in range(len(orig)):
                    prob = prob*gauss.pdf(new[i]-orig[i])
                return prob
        
        # Simple independent PDF combination
        def prior_func(priors,params):
            prob = 1.0
            for i,distr in enumerate(priors):
                prob = prob * distr.pdf(params[i])
            return prob

        f = open('ABCPredPotassium_alphabeta.txt','w')

        depols = [109, 100, 88, 76, 63, 51, 38, 32, 26, 19, 10]

        for i,d in enumerate(depols):
            print "Depolarization: "+str(d)
            [x,y,alpha,beta] = HodgkinHuxley.fig3(i)
            proto.SetInput('depolarization',V.Simple(d))
            result = fitting.approx_bayes_smc_adaptive(init,priors,x,y,prior_func,kern,distance,100,10000,0.003)

            # Write outputs to ABCPredPotassium_alphabeta.txt
            f.write(str(d)+"\n")
            f.write(str(result.pool)+"\n")
            f.write(str(result.getmean())+"\n")
            f.write(str(result.getvar())+"\n")


    # Fits Sodium conductance gating rates alpha_m, beta_m, alpha_h, and beta_h from the 
    #   simplified Hodgkin-Huxley model over a range of depolarization values.
    def TestGNaFitting(self):
        proto = fc.Protocol('projects/AidanDaly/test/protocols/hh_voltageclamp.txt')
        proto.SetOutputFolder('HodgkinHuxleyABC_alphabeta2')
        proto.SetModel('projects/AidanDaly/hodgkin_huxley.cellml', useNumba=False)
        proto.model.SetSolver(CvodeSolver())
 
        priors = [distributions.Uniform(0,10)]*4
        init = [1.0,1.0,1.0,1.0]
 
        # ABC expects this form - sets alpha/beta, runs protocol, then returns sq_err of result
        def distance(params,times,vals):
            # Need to reset model state before evaluating again
            proto.model.ResetState()
            proto.SetInput('alpha_m',V.Simple(params[0]))
            proto.SetInput('beta_m',V.Simple(params[1]))
            proto.SetInput('alpha_h',V.Simple(params[2]))
            proto.SetInput('beta_h',V.Simple(params[3]))
            proto.Run(verbose=False,writeOut=False)
            return CheckAgainst(proto.outputEnv.LookUp('time').array,proto.outputEnv.LookUp('G_Na').array,times,vals)
        
        # Kernel function is random walk with spread of 1.0
        def kern(orig,new=None):
            gauss = distributions.Normal(0.0,1.0)
            if new == None:
                new = [param + gauss.draw() for param in orig]
                return new
            else:
                prob = 1.0
                for i in range(len(orig)):
                    prob = prob*gauss.pdf(new[i]-orig[i])
                return prob
        
        # Simple independent PDF combination
        def prior_func(priors,params):
            prob = 1.0
            for i,distr in enumerate(priors):
                prob = prob * distr.pdf(params[i])
            return prob
 
        f = open('ABCPredSodium_alphabeta.txt','w')
 
        depols = [109, 100, 88, 76, 63, 51, 38, 32, 26, 19, 10, 6]
 
        for i,d in enumerate(depols):
            print "Depolarization: "+str(d)
            [x,y,alpham,betam,alphah,betah] = HodgkinHuxley.fig6(i)
            proto.SetInput('depolarization',V.Simple(d))
            print proto.inputEnv.LookUp('depolarization').unwrapped
            result = fitting.approx_bayes_smc_adaptive(init,priors,x,y,prior_func,kern,distance,100,10000,0.003)
 
            # Write outputs to ABCPredSodium_alphabeta.txt
            f.write(str(d)+"\n")
            f.write(str(result.pool)+"\n")
            f.write(str(result.getmean())+"\n")
            f.write(str(result.getvar())+"\n")


'''
    HELPER METHODS
'''

# Evaluates RMSE between experimental and predicted values
# Uses time points in simulation that are closest to experimental
def CheckAgainst(self, predTimes, predVals, experTimes, experVals):
    curr = 0
    predValsClosest = []
        
    # Finds experimental output at times closest to experimental
    for tval in experTimes:
        while tval > predTimes[curr+1]:
            if curr >= len(predTimes)-2:
                break
            curr = curr+1
        if abs(tval-predTimes[curr]) < abs(tval-predTimes[curr+1]):
            predValsClosest = predValsClosest + [predVals[curr]]
        else:
            predValsClosest = predValsClosest + [predVals[curr+1]]
    # Calculate squared error
    sq_err = 0

    for i,val in enumerate(experVals):
        sq_err = sq_err + pow(val-predValsClosest[i],2)
    return pow(sq_err/len(experVals),0.5)
    