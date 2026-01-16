#!/usr/bin/env python3

# 8 parameter scan with random values in an interval 
# Bisection method to scan over TanBeta values, only one TanBeta value per point
# Run in directory /Vevacious 



import sys # Allows to read the script name and arguments in sys.argv
import subprocess # To interact with the OS
import numpy as np 
#import math # C math standard functions



class ScanPathTo : 
 '''
 Paths to relevant files and directories.  
 '''
 MainDir = "/home/esau/Dokumente/Omar/HiggsPortals" # Contains the Vevacious/ and SPheno/ (and /Micromegas) 
 # Some files required to run the script
 LHin = MainDir + "/SPheno/Input_Files/LesHouches.in.2HSDM_low" # SPheno input
 LHspc = MainDir + "/Vevacious/SPheno.spc.2HSDM" # SPheno output 
 Vevoutxml = MainDir + "/Vevacious/results/2HSDM.vout" # Vevacious output 
 SPhenoModel = MainDir + "/SPheno/bin/SPheno2HSDM" # SPheno Model executable 
 VevPP = MainDir + "/Vevacious/bin/VevaciousPlusPlus" # Vevacious executable 
 VevInxml = MainDir + "/Vevacious/bin/2HSDMInput.xml" # Vevacious Input.xml file 
 ScanOut = "scan.dat" # Scan output file name 
 MicroDir = "/home/esau/Dokumente/micromegas_5.3.35/2HDM" # Micromegas model directory 
 Microexe = MicroDir + "/CalcOmega_MOv5" # Micromegas model executable 
 Omegaout = MainDir + "/Vevacious/omg.out" # Micromegas output 
 ScanOutMicro = "scan_stableplusrelic.dat" # Scan output file name 


class ScanVar : 
 '''
 Some global variables that define details of the parameter scan.
 '''
 vEWpdg = 246.219640 # GeV, check [PDG (2022) from Fermi coupling] v_EW == 2.46219640(63)E+02 GeV (consistent)
 #vEWerr = 0.00063 # GeV, error in the determination of v_EW
 vEWSigma = 0.000063
 tberr = 5 # Check v2/v1 == TanBeta (consistency between Vevacious results and input)
 mHiggspdg = 125.25 # Higgs mass [PDG (2022)] m_h = (125.25 \pm 0.17) GeV
 #mHiggserr = 1.7 # 5 GeV  
 mHiggsSigma = 0.17
 Hphystag = 1 
 Omgh2 = 0.120 # DM relic density [arXiv:1807.06209] Omega*h^2 = 0.120 \pm 0.001 
 #Omgh2err = 0.01
 Omgh2Sigma = 0.001 



#------------------Functions---------------------
def makeLesHouchesin(pspoint) : # pspoint is an array of floats: [Lam1,...,Lam7,Mu122in,TanBetain] 
 '''
 Writes new LesHouches.in.2HSDM_low from the original
 with given parameters (list of floats). 
 '''
 fin = open(ScanPathTo.LHin,"r") # /home/omarset/spheno/input/LesHouches.in.2HSDM_low
 dataout = ""
 replaced = False 
 count = 0

 for linea in fin : 
     for word in linea.split() : 
         if word.startswith("#") : break 
         elif word == "MINPAR" :
             # Write new Block MINPAR to LesHouches.in (numbers in scientific notation)
             dataout += "Block MINPAR      # Input parameters \n" \
             + " 1   " + str(format(pspoint[0],".7E")) + "    # Lambda1Input\n" \
             + " 2   " + str(format(pspoint[1],".7E")) + "    # Lambda2Input\n" \
             + " 3   " + str(format(pspoint[2],".7E")) + "    # Lambda3Input\n" \
             + " 4   " + str(format(pspoint[3],".7E")) + "    # Lambda4Input\n" \
             + " 5   " + str(format(pspoint[4],".7E")) + "    # Lambda5Input\n" \
             + " 6   " + str(format(pspoint[5],".7E")) + "    # Lambda6Input\n" \
             + " 7   " + str(format(pspoint[6],".7E")) + "    # Lambda7Input\n" \
             + " 9   " + str(format(pspoint[7],".7E")) + "    # M122Input\n" \
             + " 10   " + str(format(pspoint[8],".7E")) + "    # TanBeta\n" 
             # "{:.7E}".format(Lam1) equivalent
             replaced = True
             count += 1
             break
             
     if not replaced : dataout += linea
     else : 
         if count > 10 : # 10 lines are edited
             dataout += linea
             replaced = False 
         count += 1
 fin.close()
 
 fin = open(ScanPathTo.LHin,"w")
 fin.write(dataout)
 fin.close()

 del dataout 
 return 0 
 
 
 
def readSPhenospc() : 
    '''
    Reads SPheno.spc file blocks to a string, then takes only 
    relevant values to a line in a string 
    '''
    # Read SPheno output in blocks 
    fspc = open(ScanPathTo.LHspc,"r") # "/home/omarset/Vevacious/SPheno.spc.2HSDM"
    
    blocksout = ""
    blockMIN = False 
    blockMASS = False
    blockTLU = False
    blockTLUwT = False
    # Read Higgs mass and unitarity constraints 
    for linea in fspc : 
        for word in linea.split() : 
            if word.startswith("#") : break # If it is a comment go to next line 
            else :  
                if word == "MINPAR" : 
                    blockMIN = True
                    count = 1 
                if word == "MASS" : 
                    blockMASS = True
                    count = 1 
                if word == "TREELEVELUNITARITY" : 
                    blockTLU = True
                    count = 1
                if word == "TREELEVELUNITARITYwTRILINEARS" : 
                    blockTLUwT = True
                    count = 1
        if blockMIN : 
            if 1 < count < 16 : blocksout += linea 
            elif count > 15 : blockMIN = False 
            count += 1
        if blockMASS : 
            if 2 < count < 8 : blocksout += linea 
            elif count > 7: blockMASS = False 
            count += 1
        if blockTLU : 
            if count == 2 : blocksout += linea 
            elif count > 2 : blockTLU = False 
            count += 1
        if blockTLUwT : 
            if count == 2 : blocksout += linea 
            elif count > 2 : blockTLUwT = False 
            count += 1
    
    fspc.close()
    # SPheno output in a line 
    lineout = "" 
    for linea in blocksout.split("\n") :
        for index,word in enumerate(linea.split()) : 
            if word.startswith("#") : break
            elif index == 1 : lineout += word + "\t" 
    del blocksout 
    
    return lineout 
    


def readVevout() : 
    '''
    Reads the VevaciousPlusPlus output file  to a 
    string that contains a single line
    '''
    fvev = open(ScanPathTo.Vevoutxml,"r") # "/home/omarpf/HiggsPortals/Vevacious/results/2HSDM.vout"
    
    results = ""
    readnext = False
    fvals = False  
    count = 0 

    for linea in fvev : 
        if readnext : 
            results += linea.split()[0] + "\t"
            readnext = False 
        for word in linea.split() : 
            if word == "<StableOrMetastable>" : readnext = True
            if word == "<DsbVacuum>" : fvals = True 
            if word == "<DsbSurvivalProbability>" : readnext = True
            if word == "<DsbLifetime>" : readnext = True
        if fvals : 
            count += 1 
            if 2 < count < 5 : results += linea.split()[0] + "\t" 
            elif count > 4 : fvals = False 
    
    fvev.close()
    
    return results



def runSPheno() : 
 '''
 Runs SPheno, reads output to a string. 
 '''
 pspoint = ""
 print("Running SPheno...")
 
 try: 
     sphpr = subprocess.run([ScanPathTo.SPhenoModel,ScanPathTo.LHin],
         stdout = subprocess.PIPE, stderr = subprocess.PIPE, # capture_output = True,
         timeout = 4)           # python > 3.7   
 except subprocess.TimeoutExpired : 
     print("SPheno ran for too long (timeout).")
     pspoint += "SPheno exception (timeout) \t"
     return pspoint 
 
 if sphpr.returncode == 0 and sphpr.stderr == b'' :
     print("SPheno calculated the spectrum, reading output...")
     SPhenout = readSPhenospc()
     pspoint += SPhenout 
     del SPhenout 
 else : 
     print("SPheno error, code: " + str(sphpr.returncode) + ". No output.") 
     print("STDOUT: " + str(sphpr.stdout,'utf-8'))
     print("STDERR: " + str(sphpr.stderr,'utf-8'))
     pspoint += "SPheno exception, error code: " + str(sphpr.returncode) + "\t" #\
         #+ "STDOUT: " + str(sphpr.stdout).strip("b") #\
         #+ "STDERR: " + str(sphpr.stderr)
 
 return pspoint 
 
 
 
 
def runVevacious() : 
 '''
 Runs Vevacious, reads output to a string. 
 '''
 pspoint = ""
 print("Running Vevacious...")
 
 try :
     vevpr = subprocess.run([ScanPathTo.VevPP,ScanPathTo.VevInxml],
             stdout = subprocess.PIPE, stderr = subprocess.PIPE, # capture_output = True,
             timeout = 4) # python > 3.7 
 except subprocess.TimeoutExpired : 
     print("Vevacious ran for too long (timeout).")
     pspoint += "Failed to minimize (timeout) \t"
     return pspoint 
 
 if vevpr.returncode == 0 : 
     print("Vevacious found a minimum, reading output...")
     Vevout = readVevout()
     pspoint += Vevout
     del Vevout 
 elif vevpr.returncode == -6 : # forrtl: severe (174): SIGSEGV, segmentation fault occurred
     print("Vevacious failed to find a minimum (segmentation fault).")
     print(str(vevpr.stderr,'utf-8'))
     pspoint += "Failed to minimize (segmentation fault) \t" 
 else : 
     print("Vevacious failed to find a minimum (unknown error), code: "+ str(vevpr.returncode))
     print("STDOUT: " + str(vevpr.stdout,'utf-8'))
     print("STDERR: " + str(vevpr.stderr,'utf-8'))
     pspoint += "Failed to minimize (unknown error), code: " + str(vevpr.returncode) + "\t" 
 
 return pspoint 



def ViaPointC(point, fullcond = True) : 
 '''
 Check conditions for viability of parameter space point,
 use fullcond=False if there is no Vevacious output to check.
 '''
 pl = point.split() # String to a list 
 
 if pl[0] == "SPheno" : return False 
 if not fullcond : return True # Up to here checks for only SPheno output 
 
 if pl[21] == "Failed" : return False 
 unit1 = float(pl[19]) 
 unit2 = float(pl[20]) 
 if unit1 != 1 or unit2 != 1 : return False 
 # Vevacious output 
 v1 = float(pl[22]) 
 v2 = float(pl[23]) 
 vEWobs = (v1**2 + v2**2)**(1/2) 
 tbvout = v2/v1 
 tbin = float(pl[8])
 Mhh1 = float(pl[14])
 Mhh2 = float(pl[15])
 
 if pl[21] != "stable" : return False # Could add metastable points with reasonable survival probability 
 if not abs(tbin - tbvout) <= ScanVar.tberr : return False # Consistency 
 if not abs(ScanVar.vEWpdg - vEWobs) <= ScanVar.vEWerr : return False 
 if not abs(Mhh2 - Mhh1) >= 50 : return False 
 
 return True 



def calcChiVec(pspoint) : 
 '''
 Reads the parameter space point from a string, 
 calculates the chi vector. 
 '''
 pspl = pspoint.split() 
 mHiggsobs = float(pspl[14])
 v1 = float(pspl[22]) 
 v2 = float(pspl[23]) 
 vEWobs = (v1**2 + v2**2)**(1/2) 
 
 chi = [(mHiggsobs - ScanVar.mHiggspdg)/(ScanVar.mHiggsSigma), (vEWobs - ScanVar.vEWpdg)/(ScanVar.vEWSigma)]
 
 return chi 
 
 
 
def calcChiSqrd(pspoint) : 
 '''
 Reads the parameter space point from a string, 
 calculates the chi squared (Higgs mass and EW VEV), 
 adds it at the end of the string. 
 '''
 chi = calcChiVec(pspoint)
 chisqrd = chi[0]**2 + chi[1]**2 
 pspoint += str(format(chisqrd,".8E")) + "\t" 
 
 return pspoint 



def scanTanBeta(params, tbmin, tbmax) : # params is a numpy array with 8 parameters
 '''
 Bisection method applied to TanBeta parameter scanning,
 checks against Higgs mass measured value, 
 checks for v_EW and TanBeta values consistency. 
 Result to a string. 
 '''
 
 Nmaxit = 15 # Max number of iterations (they are highly time-expensive) (20)
 Nmaxstep = 10 # Max number of bisection method steps (10)
 a = tbmin; b = tbmax 
 it = 1 # Number of iterations
 step = 1 # Number of bisection method steps
 
 print("Bisection method:")
 
 while step <= Nmaxstep and it <= Nmaxit : 
  if it >= Nmaxit + 1 : break 
  
  print("Iteration " + str(it))
  
  if step == 1 : 
      params = np.append(params,a) 
      makeLesHouchesin(params) 
      pointa = runSPheno() 
      if ViaPointC(pointa,fullcond = False) : 
          mha = float(pointa.split()[14])
      else : break # Could be... else : a += 0.1; it += 1; continue 
  
      if ScanVar.mHiggspdg - ScanVar.mHiggserr <= mha <= ScanVar.mHiggspdg + ScanVar.mHiggserr : 
          return pointa 
  
  params = np.delete(params,-1)
  params = np.append(params,b) 
  makeLesHouchesin(params) 
  pointb = runSPheno()
  if ViaPointC(pointb, fullcond=False) : 
      mhb = float(pointb.split()[14]) 
  else : break # else : b -= 0.1; it += 1; continue 
  
  if ScanVar.mHiggspdg - ScanVar.mHiggserr <= mhb <= ScanVar.mHiggspdg + ScanVar.mHiggserr : 
      return pointb 
  
  tbnew = (a + b)/2
  
  params = np.delete(params,-1)
  params = np.append(params,tbnew) 
  makeLesHouchesin(params) 
  pointc = runSPheno()
  if ViaPointC(pointc, fullcond=False) : 
      mhc = float(pointc.split()[14]) 
  else : break # else : tbnew += 0.1; it += 1; continue  
  
  if ScanVar.mHiggspdg - ScanVar.mHiggserr <= mhc <= ScanVar.mHiggspdg + ScanVar.mHiggserr : 
      return pointc 
  
  step += 1; it += 1
  
  if mhc < ScanVar.mHiggspdg < mhb or mhb < ScanVar.mHiggspdg < mhc : 
      a = tbnew
  else : 
      b = tbnew
      
  if step >= Nmaxstep + 1 : break 
  
 return 0 # Bisection method failed for the point 
 
 
 
def scanpointonce(params, TanBetamin = 4, TanBetamax = 100) : 
 '''
 Evaluates the TanBeta scan for a point given by parameters 
 in a numpy array, only in one interval. Returns string. 
 '''
 print("Scanning tanbeta values " + str(TanBetamin) + " to " + str(TanBetamax) + "...")
 
 pspoint = scanTanBeta(params,TanBetamin,TanBetamax)
 if pspoint == 0 : return 0 
 if not ViaPointC(pspoint,fullcond = False) : return 0 
 
 print("Found Higgs mass value for some tanbeta.")
 
 Vevout = runVevacious()
 pspoint += Vevout 
 print("Testing consistency and viability...")
 if not ViaPointC(pspoint,fullcond = True) : return 0  
 pspoint = calcChiSqrd(pspoint)
 
 print("Found viable point!\n")
 
 return pspoint 
 
 
 
def scanpoint(params) : #, chisqrd = True) : 
 '''
 Scan for TanBeta in 2 intervals (4,100) and (100,200). 
 Returns chi squared value. 
 '''
 
 pspoint = scanpointonce(params,4,100) 
 if pspoint != 0 : return calcChiVec(pspoint) 
 
 print("First scan unsuccesful...\n") 
 
 pspoint = scanpointonce(params,100,200) 
 if pspoint != 0 : return calcChiVec(pspoint) 
 
 print("Second scan unsuccesful...\n") 
 
 return [np.float64(1E+50), np.float64(1E+50)] 
 
 
 
def modelTanBetadependence() : 
 '''
 Evaluate a specific model with SPheno for several TanBeta values. 
 -----------------  (Under construction)  ---------------------- 
 '''
 
 fout = open(ScanPathTo.ScanOut,"w") # Output file with points info. 
 
 fout.write("# Lambda1\tLambda2\tLambda3\tLambda4\tLambda5\tLambda6\tLambda7\t" 
     + "Mu122\tTanBeta\tMS2\tLambdaS\tLambda1S\tLambda2S\tLambda12S\t" 
     + "Mhh_1\tMhh_2\tMAh_2\tMHm_2\tMss\t" + "TreeLvlUnit\tTreeLvlUnitwTri\t")
 
 point = np.zeros(9)
 point[0] = 9.47409440E-03 # Lam1
 point[1] = 8.45939090E-02 # Lam2
 point[2] = 3.38621350E-01 # Lam3
 point[3] = -5.33376830E-04 # Lam4
 point[4] = 4.65478530E-04 # Lam5
 point[5] = 2.16774600E-05 # Lam6
 point[6] = 7.61827300E-05 # Lam7 
 point[7] = -8.37928860E+04 # Mu122
 
 Npoints = 100
 nviable = 0
 for n in range(Npoints + 1) :
  if n < 80 : point[8] = 1 + (n/80)*(50-1) # point[8] = TanBeta
  else : point[8] = 50 + ((n-80)/21)*(200-50)
  makeLesHouchesin(point)
  pointout = runSPheno()
  if ViaPointC(pointout, fullcond = False) :
      nviable += 1
      fout.write(pointout + "\n")
 fout.close()
 
 print("\nFinished parameter scan!\n")
 print("Tested " + str(Npoints) + " random parameter space points.")
 print("Found " + str(nviable) + " viable points.")
 print("Wrote results to file " + ScanPathTo.ScanOut + ".\n")
 
 return 0 
 
 
 
def randscan() : 
 '''
 Scan with random parameter input in a defined 
 interval read from input file. 
 '''
 print("Parameter Scan (v. 1.1)\n")
 
 datlist = np.genfromtxt(sys.argv[1], usecols = (1), max_rows = 17) # Read input file 
 # Asign the variables
 Lam1min = datlist[0]
 Lam1max = datlist[1]
 Lam2min = datlist[2]
 Lam2max = datlist[3]
 Lam3min = datlist[4]
 Lam3max = datlist[5]
 Lam4min = datlist[6]
 Lam4max = datlist[7]
 Lam5min = datlist[8]
 Lam5max = datlist[9]
 Lam6min = datlist[10]
 Lam6max = datlist[11]
 Lam7min = datlist[12]
 Lam7max = datlist[13]
 Mu122min = datlist[14]
 Mu122max = datlist[15]
 Npoints = int(datlist[16])
 del datlist # Delete list
 
 fout = open(ScanPathTo.ScanOut,"w") # Output file with points info. 

 fout.write("# Lambda1\tLambda2\tLambda3\tLambda4\tLambda5\tLambda6\tLambda7\t" 
     + "Mu122\tTanBeta\tMS2\tLambdaS\tLambda1S\tLambda2S\tLambda12S\t" 
     + "Mhh_1\tMhh_2\tMAh_2\tMHm_2\tMss\t" + "TreeLvlUnit\tTreeLvlUnitwTri\t"
     + "Stability\tFieldValues: v_1\tv_2\tChi squared\n")
     #Survival Probability\tLifetime\n")

 nviable = 0 
 rng = np.random.default_rng() 
 
 for n in range(Npoints) :
  params = np.zeros(8)
  # Logarithmic steps 
  params[0] = Lam1min*(Lam1max/Lam1min)**(rng.random()) # Lam1 
  params[1] = Lam2min*(Lam2max/Lam2min)**(rng.random()) # Lam2
  params[2] = Lam3min*(Lam3max/Lam3min)**(rng.random()) # Lam3
  params[3] = Lam4min*(Lam4max/Lam4min)**(rng.random()) # Lam4
  params[4] = Lam5min*(Lam5max/Lam5min)**(rng.random()) # Lam5
  params[5] = Lam6min*(Lam6max/Lam6min)**(rng.random()) # Lam6
  params[6] = Lam7min*(Lam7max/Lam7min)**(rng.random()) # Lam7 
  params[7] = Mu122min*(Mu122max/Mu122min)**(rng.random()) # Mu122
  #params[5] = 0 # Lam6 = 0
  #params[6] = 0 # Lam7 = 0
  #params[7] = 0 # Mu122 = 0
  
  print("Parameter space point " + str(n + 1) + "/" + str(Npoints) + "\n")
  
  if (n+1)%10 == 0 : fout.flush() # Empty RAM to disk 
  
  pspoint = scanpointonce(params,4,100) 
  if pspoint != 0 : 
      fout.write(pspoint + "\n")
      nviable += 1 
      print("Added to results. ("+ str(nviable) + "/" + str(n + 1) + ")\n")
      onemorepls = False 
  onemorepls = True 
  
  if not onemorepls : continue 
  
  print("First scan unsuccesful...\n") 
 
  pspoint = scanpointonce(params,100,200) 
  if pspoint != 0 : 
      fout.write(pspoint + "\n")
      nviable += 1 
      print("Added to results. ("+ str(nviable) + "/" + str(n + 1) + ")\n")
  
  print("Second scan unsuccesful...\n") 
  
  print("Non-viable point. (" + str(nviable) + "/" + str(n + 1) + " viable points so far)\n")
  
 fout.close() 
 
 print("\nFinished parameter scan!\n")
 print("Tested " + str(Npoints) + " random parameter space points.")
 print("Found " + str(nviable) + " viable points.")
 print("Wrote results to file " + ScanPathTo.ScanOut + ".\n")
 
 return 0 
 
 
 
def main() : 
 '''
 Main funciton. 
 ''' 
 
 #modelTanBetadependence()
 randscan() 
 
 return 0 



if __name__ == '__main__' : 
    main()
    
    
    
