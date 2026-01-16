#!/usr/bin/env python3

# Version 0.1 
# Chi squared minimizer 
# Run in directory /Vevacious 



import numpy as np 
import time 
import lmfit # https://lmfit.github.io/lmfit-py/
from lmfit import minimize, Parameters, fit_report
import emcee # MCMC
#from iminuit import Minuit 

import spvevmicro as svm 


class MinimizeStopper : 
 '''
 Behaves as a function to call after each objective function iteration. 
 '''
 def __init__(self, max_sec = 60) : # Runs when instance is created 
     self.max_sec = max_sec 
     #self.start = time.time() #time.time_ns() # python3.7 and older 
 def __call__(self, params, iter, resid) : # Runs when is called (as function) 
     print(f"Iteration {iter}") # current iteration number 
     print(params.valuesdict()) # current parameters 
     print(f"Current residual: {resid}") # current chi vector
     chisqrd = np.inner(resid,resid)
     print(f"Chi squared: {chisqrd}")
     #elapsed = time.time() - self.start 
     #if elapsed > self.max_sec : 
         #return True 
 
 
 
def scan_lmfit_p1(max_time = 60*60*2, max_fncev = 7200) : # 1 hour timeout ~ 3600 iterations
 '''
 Paso 1. 
 Scan 14 parameters. 
 Using 10sigma intervals: vEWerr = 0.00063, mHiggserr = 1.7, Omgh2err = 0.01 and tberr = 1. 
 '''
 '''
 pspoint_diffevo = Parameters()
 # Best point old scan (chi^2 = 40.53)
 pspoint_diffevo.add('Lam1', value = 0.00026921048, min = 1E-04, max = 4E-04, vary = False)
 pspoint_diffevo.add('Lam2', value = 0.12328174, min = 0.1, max = 0.14, vary = False)
 pspoint_diffevo.add('Lam3', value = 0.49571412, min = 0.2, max = 0.5, vary = False)
 pspoint_diffevo.add('negLam4', value = 0.15779285, min = 0.1, max = 0.2, vary = False)
 pspoint_diffevo.add('Lam5', value = 0.0079333098, min = 5E-03, max = 1E-02, vary = False)
 pspoint_diffevo.add('Lam6', value = 0.00042895379, min = 1E-08, max = 1E-03, vary = False)
 pspoint_diffevo.add('Lam7', value = 0.0002791794, min = 1E-08, max = 1E-03, vary = False)
 pspoint_diffevo.add('negMu122', value = 1462.8741151948557, min = 1E+02, max = 1E+04, vary = True)
 pspoint_diffevo.add('TanBeta', value = 102.651136073714, min = 4, max = 200, vary = False)
 pspoint_diffevo.add('MS2', value = 78325.53039861172, min = 8E+03, max = 1E+06, vary = True) 
 pspoint_diffevo.add('LamS', value = 0.006142983436304768, min = 0.005, max = 0.15, vary = False)
 pspoint_diffevo.add('Lam1S', value = 0.49631019904995227, min = 0.0001, max = 0.5, vary = False)
 pspoint_diffevo.add('Lam2S', value = 0.11836278048716235, min = 0.0001, max = 0.4, vary = False)
 pspoint_diffevo.add('Lam12S', value = 0.002444277324586075, min = 0.000001, max = 0.1, vary = False)
 '''
 pspoint_leastsquares = Parameters() 
 # Chi^2 ~ 10^3
 pspoint_leastsquares.add('Lam1', value = 0.00026921048, min = 1E-04, max = 4E-04) #, vary = False)
 pspoint_leastsquares.add('Lam2', value = 0.12328174, min = 0.1, max = 0.14)
 pspoint_leastsquares.add('Lam3', value = 0.49571412, min = 0.1, max = 0.6)
 pspoint_leastsquares.add('negLam4', value = 0.15779285, min = 0.01, max = 0.4)
 pspoint_leastsquares.add('Lam5', value = 0.0079333098, min = 5E-03, max = 1E-02)
 pspoint_leastsquares.add('Lam6', value = 0.00042895379, min = 1E-04, max = 1E-03)
 pspoint_leastsquares.add('Lam7', value = 0.0002791794, min = 8E-05, max = 1E-03)
 pspoint_leastsquares.add('negMu122', value = 1462.8741151948557, min = 1.0E+02, max = 1.0E+04)
 pspoint_leastsquares.add('TanBeta', value = 102.651136073714, min = 4, max = 200)
 pspoint_leastsquares.add('MS2', value = 78325.53039861172, min = 1E+03, max = 1E+06)
 pspoint_leastsquares.add('LamS', value = 0.006142983436304768, min = 0.0001, max = 0.1)
 pspoint_leastsquares.add('Lam1S', value = 0.49631019904995227, min = 0.0001, max = 0.5)
 pspoint_leastsquares.add('Lam2S', value = 0.11836278048716235, min = 0.0001, max = 0.3)
 pspoint_leastsquares.add('Lam12S', value = 0.002444277324586075, min = 0.00001, max = 0.1)
 
 #print(pspoint.valuesdict()['Lam1']) es 2.69210480E-04
 
 result = minimize(svm.runSVM, params = pspoint_leastsquares, 
         method = 'least_squares', # 'differential_evolution' (uses the interval, not the initial point)
         iter_cb = MinimizeStopper(max_time), # Max time in seconds 
         max_nfev = max_fncev) # Max number of calls to objective function 
         #ftool=1e-3, xtool=1e-3) # **fit_kws
 
 print(f"Success: {result.success}")
 if result.success :  
     #print(f"Result status: {result.status}")
     print(f"Message: {result.message}")
 #result.params.pretty_print()
 print(fit_report(result))
 
 return 0 
 
 
def scan_lmfit_p2(max_time = 60*60*2) : 
 '''
 Paso 2. 
 '''
 
 point = Parameters() 
 # Chi^2 == 10.63 point with real good characteristics 
 point.add('Lam1', value = 0.00026921048, min = 1E-04, max = 4E-04, vary = True)  
 point.add('Lam2', value = 0.12328174, min = 0.1, max = 0.14, vary = True)
 point.add('Lam3', value = 0.49571412, min = 0.1, max = 0.6, vary = True)
 point.add('negLam4', value = 0.15779285, min = 0.01, max = 0.4, vary = True)
 point.add('Lam5', value = 0.0079333098, min = 5E-03, max = 1E-02, vary = True)
 point.add('Lam6', value = 0.00042895379, min = 1E-04, max = 1E-03, vary = True)
 point.add('Lam7', value = 0.0002791794, min = 8E-05, max = 1E-03, vary = True)
 point.add('negMu122', value = 1462.8741151948557, min = 1.0E+02, max = 1.0E+04, vary = True)
 point.add('TanBeta', value = 102.651136073714, min = 4, max = 200, vary = True)
 point.add('MS2', value = 78325.53039861172, min = 1E+03, max = 1E+06, vary = True)
 point.add('LamS', value = 0.006142983436304768, min = 0.0001, max = 0.1, vary = True)
 point.add('Lam1S', value = 0.49631019904995227, min = 0.0001, max = 0.5, vary = True)
 point.add('Lam2S', value = 0.11836278048716235, min = 0.0001, max = 0.3, vary = True)
 point.add('Lam12S', value = 0.002444277324586075, min = 0.00001, max = 0.1, vary = True)

 result = minimize(svm.runSVM, params = point, 
         method = 'emcee', 
         nan_policy = 'omit', steps = 1000, #thin = 5, burn = 0, 
         nwalkers = 160, 
         is_weighted = True,
         iter_cb = MinimizeStopper(max_time)) 
 
 print(f"Success: {result.success}")
 if result.success :  
     print(f"Message: {result.message}")
     
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
 
 params_omar1 = { # chi^2 ~ 2000, tbvout == 102.080094199
     'Lam1': 0.00028439245, 
     'Lam2': 0.11592538801699866, 
     'Lam3': 0.48687466, 
     'negLam4': 0.14046414, 
     'Lam5': 0.0083945649, 
     'Lam6': 0.00048051592, 
     'Lam7': 0.00029747456, 
     'negMu122': 1613.8836, 
     'TanBeta': 101.77094, 
     'MS2': 403.08351, 
     'LamS': 0.11481724, 
     'Lam1S': 0.11766687, 
     'Lam2S': 0.10692459976849451, 
     'Lam12S': 0.14860111
 }
 
 params_esau2 = { # Chi^2 = 9.2642, abs(tbin-tbvout) = 0.098657494
# 2.69210480E-04	1.23281740E-01	4.95714120E-01	-1.57792850E-01	7.93330980E-03	4.28953790E-04	2.79179400E-04	-1.58424950E+03	1.02651140E+02	2.48429900E+02	7.01381850E-02	2.17325910E-01	1.36131800E-01	9.48497540E-04	1.24814657E+02	3.98817527E+02	3.97526817E+02	4.03331152E+02	6.83603171E+01	1.00000000E+00	1.00000000E+00	stable	2.4008	246.208	0.118718	9.26422622E+00	
# Current residual: [-2.5608411764706065, 1.0309193399619303, -1.2819999999999914]
     'Lam1': 0.00026921048, 
     'Lam2': 0.12328174, 
     'Lam3': 0.49571412, 
     'negLam4': 0.15779285, 
     'Lam5': 0.0079333098, 
     'Lam6': 0.00042895379, 
     'Lam7': 0.0002791794, 
     'negMu122': 1584.2495140779977, 
     'TanBeta': 102.651136073714, 
     'MS2': 248.42990354079654, 
     'LamS': 0.07013818516136634, 
     'Lam1S': 0.21732590735076793, 
     'Lam2S': 0.13613179724942037, 
     'Lam12S': 0.00094849754198518
 }
 
 params_esau3 = { # Chi squared = 7.66, abs(tbin-tbvout) = 0.09438572
# 2.69210480E-04	1.23281740E-01	4.95714120E-01	-1.57792850E-01	7.93330980E-03	4.28953790E-04	2.79179400E-04	-1.58362160E+03	1.02651140E+02	1.18669250E+02	6.14298340E-03	3.88703910E-01	1.29583430E-01	9.94951310E-04	1.24813227E+02	3.98734352E+02	3.97443391E+02	4.03249886E+02	6.79129685E+01	1.00000000E+00	1.00000000E+00	stable	2.4007	246.208	0.120186	7.66678005E+00	
# Current residual: [-2.569252941176484, 1.0154424537342983, 0.1860000000000056]
     'Lam1': 0.00026921048, 
     'Lam2': 0.12328174, 
     'Lam3': 0.49571412, 
     'negLam4': 0.15779285, 
     'Lam5': 0.0079333098, 
     'Lam6': 0.00042895379, 
     'Lam7': 0.0002791794, 
     'negMu122': 1583.6216110734883, 
     'TanBeta': 102.651136073714, 
     'MS2': 9448.658514581215, 
     'LamS': 0.006142983436304768, 
     'Lam1S': 0.3887039091005354, 
     'Lam2S': 0.12958342610756332, 
     'Lam12S': 0.000994951314220159
 }
 
 params_omar2 = { # Chi squared = 10.63, abs(tbin-tbvout) = 0.466824077
# 2.69210480E-04	1.23281740E-01	4.95714120E-01	-1.57792850E-01	7.93330980E-03	4.28953790E-04	2.79179400E-04	-1.46287410E+03	1.02651140E+02	7.83255300E+04	6.14298340E-03	4.96310200E-01	1.18362780E-01	2.44427730E-03	1.24870023E+02	3.82344328E+02	3.81004591E+02	3.87346524E+02	2.87253596E+02	1.00000000E+00	1.00000000E+00	stable	2.40945	246.208	0.119904	1.06320540E+01
# Current residual: [-2.235158823529392, 2.3721094102877522, -0.09599999999999886]
     'Lam1': 0.00026921048, 
     'Lam2': 0.12328174, 
     'Lam3': 0.49571412, 
     'negLam4': 0.15779285, 
     'Lam5': 0.0079333098, 
     'Lam6': 0.00042895379, 
     'Lam7': 0.0002791794, 
     'negMu122': 1462.8741151948557,
     'TanBeta': 102.651136073714, 
     'MS2': 78325.53039861172, 
     'LamS': 0.006142983436304768, 
     'Lam1S': 0.49631019904995227, 
     'Lam2S': 0.11836278048716235, 
     'Lam12S': 0.002444277324586075
 }
 
 result = svm.runSVM(params_omar2, False)
 
 print(result[0])
 print(f"Current residual: {result[1]}") #print(result[1]) 
 
 return 0 



def main() : 
 '''
 Main function, choose what to do. 
 '''
 #ftest() # To use comment line 55 in spvevmicro.py 
 #scan_lmfit_p1()
 scan_lmfit_p2()
 
 return 0  
 
 
 
if __name__ == '__main__' : 
    main()
    
    
    
