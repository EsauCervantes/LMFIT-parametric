#!/usr/bin/env python3

# Version 0.1 
# Chi squared minimizer 
# Run in directory /Vevacious 



import numpy as np 
import time 
import lmfit # https://lmfit.github.io/lmfit-py/
from lmfit import minimize, Parameters, fit_report
#import emcee # MCMC
#from iminuit import Minuit 

import spvevmicro as svm 


class MinimizeStopper : 
 '''
 Behaves as a function to call after each objective function iteration. 
 '''
 def __init__(self, max_sec = 60) : # Runs when instance is created 
     self.max_sec = max_sec 
     self.start = time.time() #time.time_ns() # python3.7 and older 
 def __call__(self, params, iter, resid) : # Runs when is called (as function) 
     print(f"Iteration {iter}") # current iteration number 
     print(params.valuesdict()) # current parameters 
     print(f"Current residual: {resid}") # current chi vector
     elapsed = time.time() - self.start 
     if elapsed > self.max_sec : 
         return True 
 
 
 
def scan_lmfit_p1(max_time = 60*5, max_fncev = 250) : 
 '''
 Paso 1. 
 Scan 14 parameters. 
 Using 10sigma intervals: vEWerr = 0.00063, mHiggserr = 1.7, Omgh2err = 0.01 and tberr = 1. 
 '''
 
 pspoint = Parameters()
 # Best point old scan (chi^2 = 40.53)
 pspoint.add('Lam1', value = np.log10(2.69210480E-04), min = -8.0, max = np.log10(0.75)) #, vary = False)
 pspoint.add('Lam2', value = np.log10(1.23281740E-01), min = -1.0, max = np.log10(0.14))
 pspoint.add('Lam3', value = np.log10(4.95714120E-01), min = -8.0, max = np.log10(0.75))
 pspoint.add('negLam4', value = np.log10(1.57792850E-01), min = -8.0, max = np.log10(0.75))
 pspoint.add('Lam5', value = np.log10(7.93330980E-03), min = -8.0, max = np.log10(0.75))
 pspoint.add('Lam6', value = np.log10(4.28953790E-04), min = -8.0, max = np.log10(0.75))
 pspoint.add('Lam7', value = np.log10(2.79179400E-04), min = -8.0, max = np.log10(0.75))
 pspoint.add('negMu122', value = np.log10(1.62386970E+03), min = 0.0, max = 4.0)
 pspoint.add('TanBeta', value = np.log10(1.00000000E+02), min = np.log10(4), max = np.log10(200))
 pspoint.add('MS2', value = np.log10(5.00000000E+02), min = 2.0, max = 6.0)
 pspoint.add('LamS', value = np.log10(1.00000000E-01), min = -8.0, max = np.log10(0.75))
 pspoint.add('Lam1S', value = np.log10(1.00000000E-01), min = -8.0, max = np.log10(0.2))
 pspoint.add('Lam2S', value = np.log10(1.00000000E-01), min = -8.0, max = np.log10(0.2))
 pspoint.add('Lam12S', value = np.log10(1.00000000E-01), min = -8.0, max = np.log10(0.2))
 
 #print(pspoint.valuesdict()['Lam1']) es 2.69210480E-04
 
 result = minimize(svm.runSVM, params = pspoint, 
         method = 'least_squares', # 'differential_evolution' (big step, needs better error handling) 
         iter_cb = MinimizeStopper(max_time), # Max time in seconds 
         max_nfev = max_fncev) # Max number of calls to objective function 
         #ftool=1e-3, xtool=1e-3) # **fit_kws
 
 print(f"Success: {result.success}")
 if result.success :  
     print(f"Result status: {result.status}")
     print("Message: " + result.message)
 #result.params.pretty_print()
 print(fit_report(result))
 
 return 0 
 
 
def scan_lmfit_p2() : 
 '''
 Paso 2. 
 '''

 result = minimize(svm.runSVM, params = params, method = 'emcee', 
         nan_policy = 'omit', burn = burn, steps = steps, thin = thin, 
         nwalkers = nwalkers, is_weighted = True) 
 
 print(fit_report(result))
 
 return 0 
 
 

def ftest() : 
 '''
 Test. 
 Use: vEWerr = 10, tberr = 5,  mHiggserr = 10, Omgh2err = 1 
 '''
 #params = np.zeros(14)
 params = {
     'Lam1':2.77739690E-04, # 2.69210480E-04
     'Lam2':1.24166280E-01, # 1.23281740E-01
     'Lam3':4.97115780E-01, # 4.95714120E-01
     'negLam4':-1.61528750E-01, # -1.57792850E-01
     'Lam5':7.88514650E-03, # 7.93330980E-03
     'Lam6':4.09505660E-04, # 4.28953790E-04
     'Lam7':2.54813760E-04, # 2.79179400E-04
     'negMu122':-1.62229490E+03, # -1.62386970E+03
     'TanBeta':1.0864900E+02, # 1.00000000E+02
     'MS2':1.8500000E+05, # 5.00000000E+02
     'LamS':1.0000000E-01, # 1.00000000E-01
     'Lam1S':1.0000000E-01, # 1.00000000E-01
     'Lam2S':1.0000000E-01, # 1.00000000E-01
     'Lam12S':1.0000000E-04 # 1.00000000E-01
 }
 
 params1 = {
     'Lam1':2.69210480E-04,
     'Lam2':1.23281740E-01,
     'Lam3':4.95714120E-01,
     'negLam4':-1.57792850E-01,
     'Lam5':7.93330980E-03,
     'Lam6':4.28953790E-04,
     'Lam7':2.79179400E-04,
     'negMu122':-1.62386970E+03,
     'TanBeta':1.00000000E+02,
     'MS2':5.00000000E+02,
     'LamS':1.00000000E-01,
     'Lam1S':1.00000000E-01,
     'Lam2S':1.00000000E-01,
     'Lam12S':1.00000000E-01
 }
 
 result = svm.runSVM(params1, False)
 
 print(result[0])
 print(result[1]) 
 
 return 0 



def main() : 
 '''
 '''
 #ftest() # To use comment line 55 in spvevmicro.py 
 scan_lmfit_p1()
 return 0  
 
 
 
if __name__ == '__main__' : 
    main()
    
    
    
