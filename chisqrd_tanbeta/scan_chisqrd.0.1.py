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

import spvevmicro_tanbeta as svm 


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
 
 
 
def scan_lmfit_p1(max_time = 60*60*1, max_fncev = 2500) : # 1 hour timeout ~ 3600 iterations
 '''
 Paso 1. 
 Scan 14 parameters. 
 Using 10sigma intervals: vEWerr = 0.00063, mHiggserr = 1.7, Omgh2err = 0.01 and tberr = 1. 
 '''
 '''
 pspoint_diffevo = Parameters()
 # Best point old scan (chi^2 = 40.53)
 pspoint_diffevo.add('Lam1', value = 2.69210480E-04, min = 2E-04, max = 3E-04) #, vary = False)
 pspoint_diffevo.add('Lam2', value = 1.23281740E-01, min = 0.1, max = 0.14)
 pspoint_diffevo.add('Lam3', value = 4.95714120E-01, min = 0.48, max = 0.5)
 pspoint_diffevo.add('negLam4', value = 1.57792850E-01, min = 0.14, max = 0.16)
 pspoint_diffevo.add('Lam5', value = 7.93330980E-03, min = 7.5E-03, max = 8.5E-03)
 pspoint_diffevo.add('Lam6', value = 4.28953790E-04, min = 4E-04, max = 5E-04)
 pspoint_diffevo.add('Lam7', value = 2.79179400E-04, min = 2E-04, max = 3E-04)
 pspoint_diffevo.add('negMu122', value = 1.62386970E+03, min = 1.6E+03, max = 1.65E+03)
 pspoint_diffevo.add('TanBeta', value = 1.00000000E+02, min = 4, max = 200)
 pspoint_diffevo.add('MS2', value = 5.00000000E+02, min = 4E+02, max = 6E+02) # Maybe bigger intervals for these 
 pspoint_diffevo.add('LamS', value = 1.00000000E-01, min = 0.05, max = 0.15)
 pspoint_diffevo.add('Lam1S', value = 1.00000000E-01, min = 0.05, max = 0.15)
 pspoint_diffevo.add('Lam2S', value = 1.00000000E-01, min = 0.05, max = 0.15)
 pspoint_diffevo.add('Lam12S', value = 1.00000000E-01, min = 0.05, max = 0.15)
 '''
 pspoint_leastsquares = Parameters() 
 # Chi^2 ~ 10^3
 pspoint_leastsquares.add('Lam1', value = 2.8439245E-04, min = 2E-04, max = 4E-04) #, vary = False)
 pspoint_leastsquares.add('Lam2', value = 1.1589472E-01, min = 0.1, max = 0.14)
 pspoint_leastsquares.add('Lam3', value = 4.8687466E-01, min = 0.3, max = 0.6)
 pspoint_leastsquares.add('negLam4', value = 1.4046414E-01, min = 0.08, max = 0.2)
 pspoint_leastsquares.add('Lam5', value = 8.3945649E-03, min = 7.5E-03, max = 9.5E-03)
 pspoint_leastsquares.add('Lam6', value = 4.8051592E-04, min = 3E-04, max = 6E-04)
 pspoint_leastsquares.add('Lam7', value = 2.9747456E-04, min = 1E-04, max = 5E-04)
 pspoint_leastsquares.add('negMu122', value = 1.6138836E+03, min = 1.0E+02, max = 1.0E+04)
 pspoint_leastsquares.add('TanBeta', value = 1.0177094E+02, min = 4, max = 200)
 pspoint_leastsquares.add('MS2', value = 4.0308351E+02, min = 1E+02, max = 1E+03)
 pspoint_leastsquares.add('LamS', value = 1.1481724E-01, min = 0.01, max = 0.2)
 pspoint_leastsquares.add('Lam1S', value = 1.1766687E-01, min = 0.01, max = 0.2)
 pspoint_leastsquares.add('Lam2S', value = 1.0529379E-01, min = 0.01, max = 0.2)
 pspoint_leastsquares.add('Lam12S', value = 1.4860111E-01, min = 0.01, max = 0.2)
 #print(pspoint.valuesdict()['Lam1']) es 2.69210480E-04
 
 result = minimize(svm.runSVM, params = pspointleastsquares, 
         method = 'least_squares', # 'differential_evolution' (uses the interval, not the initial point)
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
 
 params_esau1 = { # Chi^2 ~ 10^3
     'Lam1': 0.0002843924462262504, 
     'Lam2': 0.11589471535857743, 
     'Lam3': 0.48687466436386767, 
     'negLam4': 0.1404641396742865, 
     'Lam5': 0.0083945649381763, 
     'Lam6': 0.00048051592065202113, 
     'Lam7': 0.0002974745631670949, 
     'negMu122': 1613.8835761317557, 
     'TanBeta': 101.77093898124417, 
     'MS2': 403.0835132966783, 
     'LamS': 0.11481724488377305, 
     'Lam1S': 0.11766686769452614, 
     'Lam2S': 0.10529378904461634, 
     'Lam12S': 0.14860111071917753
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
    
    
    
