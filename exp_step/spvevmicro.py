#!/usr/bin/env python3

# Version 0.1 
# Module for chi^2 minimizer scan 



import sys # Allows to read the script name and arguments in sys.argv
import subprocess # To interact with the OS
import numpy as np 


#--------------------Classes-----------------------
class ScanPathTo : 
 '''
 Paths to relevant files and directories.  
 '''
 MainDir = "/home/omar/HiggsPortals" # Contains the Vevacious/ and SPheno/ (and /Micromegas) 
 # Some files required to run the script
 LHin = MainDir + "/SPheno/input/LesHouches.in.2HSDM_low" # SPheno input
 LHspc = MainDir + "/Vevacious/SPheno.spc.2HSDM" # SPheno output 
 Vevoutxml = MainDir + "/Vevacious/results/2HSDM.vout" # Vevacious output 
 SPhenoModel = MainDir + "/SPheno/bin/SPheno2HSDM" # SPheno Model executable 
 VevPP = MainDir + "/Vevacious/bin/VevaciousPlusPlus" # Vevacious executable 
 VevInxml = MainDir + "/Vevacious/bin/2HSDMInput.xml" # Vevacious Input.xml file 
 MicroDir = MainDir + "/Micromegas/2HSDM" # Micromegas model directory 
 Microexe = MicroDir + "/CalcOmega_MOv5" # Micromegas model executable 
 Omegaout = MainDir + "/Vevacious/omg.out" # Micromegas output 
 #ScanOut = "scan.dat" # Scan output file name 
 #ScanOutMicro = "scan_stableplusrelic.dat" # Scan output file name 


class ScanVar : 
 '''
 Some global variables that define details of the parameter scan.
 '''
 vEWpdg = 246.219640; vEWSigma = 0.000063 # GeV [PDG (2022) from Fermi coupling], v_EW == 2.46219640(63)E+02 GeV 
 vEWerr = 0.00063 # GeV error in the determination of v_EW, 10sigma -> 0.00063
 tberr = 5 # Check v2/v1 == TanBeta (consistency between Vevacious results and input)
 mHiggspdg = 125.25; mHiggsSigma = 0.17 # GeV Higgs mass [PDG (2022)] m_h = (125.25 \pm 0.17) GeV
 mHiggserr = 1.7 # GeV 10sigma -> 1.7
 Hphystag = 1 # 2 then scanTanBeta() looks for light exotic neutral Higgs
 Omgh2 = 0.120; Omgh2Sigma = 0.001  # DM relic density [arXiv:1807.06209] Omega*h^2 = 0.120 \pm 0.001 
 Omgh2err = 0.01 # 10sigma -> 0.01



#------------------Functions---------------------
def makeLesHouchesin(point, exp_step = True) : # point is an lmfit.parameter.Parameter object with 14 parameters 
 '''
 Writes new LesHouches.in.2HSDM_low from the original
 with given parameters (list of floats). 
 Sets Lam1S=Lam2S and Lam12S=0 if Lam6=0=Lam7. 
 '''
 pspoint = point.valuesdict() # Dictionary of parameters 
 if exp_step : 
     for key,value in pspoint.items() : 
         pspoint[key] = np.power(10, value) 
 
 fin = open(ScanPathTo.LHin,"r") 
 dataout = ""
 replaced = False 
 count = 0

 for linea in fin : 
     for index,word in enumerate(linea.split()) :
         if word.startswith("#") : break
         elif word == "MINPAR" :
             # Write new Block MINPAR to LesHouches.in (numbers in scientific notation)
             dataout += "Block MINPAR      # Input parameters \n" \
             + " 1   " + str(format(pspoint['Lam1'],".7E")) + "    # Lambda1Input\n" \
             + " 2   " + str(format(pspoint['Lam2'],".7E")) + "    # Lambda2Input\n" \
             + " 3   " + str(format(pspoint['Lam3'],".7E")) + "    # Lambda3Input\n" \
             + " 4   -" + str(format(pspoint['negLam4'],".7E")) + "    # Lambda4Input\n" \
             + " 5   " + str(format(pspoint['Lam5'],".7E")) + "    # Lambda5Input\n" \
             + " 6   " + str(format(pspoint['Lam6'],".7E")) + "    # Lambda6Input\n" \
             + " 7   " + str(format(pspoint['Lam7'],".7E")) + "    # Lambda7Input\n" \
             + " 9   -" + str(format(pspoint['negMu122'],".7E")) + "    # M122Input\n" \
             + " 10   " + str(format(pspoint['TanBeta'],".7E")) + "    # TanBeta\n" \
             + " 11   " + str(format(pspoint['MS2'],".7E")) + "    # MS2Input\n" \
             + " 12   " + str(format(pspoint['LamS'],".7E")) + "    # LambdaSInput\n" \
             + " 13   " + str(format(pspoint['Lam1S'],".7E")) + "    # Lambda1SInput\n" \
             + " 14   " + str(format(pspoint['Lam2S'],".7E")) + "    # Lambda2SInput\n" \
             + " 15   " + str(format(pspoint['Lam12S'],".7E")) + "    # Lambda12SInput\n"
             '''
             if pspoint[5] == 0 and pspoint[6] == 0 : 
                 dataout += " 15   " + "0.0000000E+00" + "    # Lambda12SInput\n" 
             else : 
                 dataout += " 15   " + "1.0000000E-05" + "    # Lambda12SInput\n" 
             '''
             replaced = True 
             count += 1
             break 
             
     if not replaced : dataout += linea
     else : 
         if count > 15 : # 15 lines are edited
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



def runMicrOmegas() : 
 '''
 Runs MicrOmegas and returns a string which is 
 the resulting relic density read from the file omg.out. 
 '''
 subprocess.run(["cp", ScanPathTo.LHspc, ScanPathTo.MicroDir]) # Copy SPheno output to MicrOmegas directory 
 omgout = ""
 print("Running Micromegas...") 
 
 try : 
     migaspr = subprocess.run([ScanPathTo.Microexe], 
         stdout = subprocess.PIPE, stderr = subprocess.PIPE, # capture_output = True,
         timeout = 4) #  python > 3.7
 except subprocess.TimeoutExpired : 
     print("Micromegas ran for too long (timeout).")
     print("STDOUT: " + str(migaspr.stdrout,'utf-8'))
     print("STDERR: " + str(migaspr.stderr,'utf-8'))
     omgout += "Micromegas took too long (timeout)\t"
     return omgout
 
 if migaspr.returncode != 0 : 
     print("Micromegas exception, error code: " + str(migaspr.returncode))
     print("STDOUT: " + str(migaspr.stdrout,'utf-8'))
     print("STDERR: " + str(migaspr.stderr,'utf-8'))
     omgout += "Micromegas error, code: " + str(migaspr.returncode)
     return omgout
 
 print("Micromegas calculated the relic density, reading output...")
     
 fomg = open(ScanPathTo.Omegaout,"r") # Omegaout = "/home/omarset/Vevacious/omg.out"
 capture = False 
 done = False 
 # Read output 
 for line in fomg : 
  if done : break 
  for index,word in enumerate(line.split()) : 
   if word.startswith("#") : break
   if capture and index == 1 : 
       omgout += word + "\t"
       done = True 
       capture = False 
       break
   if word == "1" and index == 0 : 
       capture = True 
       continue 
    
 fomg.close()
 return omgout 



def checkSPout(point) : 
 '''
 Check SPheno output. 
 '''
 pl = point.split() # String to a list 
 if pl[0] == "SPheno" : return False # No further calculations 
 
 unit1 = float(pl[19]) 
 unit2 = float(pl[20]) 
 if unit1 != 1 or unit2 != 1 : return False 
 '''
 Mhh1 = float(pl[14]) 
 Mhh2 = float(pl[15]) 
 if ScanVar.Hphystag == 1 : 
     if not abs(Mhh2 - Mhh1) >= 50 : return False 
 elif ScanVar.Hphystag == 2 : 
     if not 85 <= Mhh1 <= 105 : return False 
 else : return False 
 ''' 
 return True # Proceed with calculations 



def checkVEVout(point) : 
 '''
 Check Vevacious output. 
 '''
 pl = point.split() 
 if pl[21] == "Failed" : return False 
 
 if pl[21] != "stable" : return False # Could add metastable points with reasonable survival probability 
 
 tbin = float(pl[8])
 v1 = float(pl[22]) 
 v2 = float(pl[23]) 
 tbvout = v2/v1 
 if not abs(tbin - tbvout) <= ScanVar.tberr : return False # Consistency 
 '''
 vEWobs = (v1**2 + v2**2)**(1/2) 
 if not abs(ScanVar.vEWpdg - vEWobs) <= ScanVar.vEWerr : return False 
 '''
 return True 
 
 
 
def checkMicrout(point) : 
 '''
 Checks relic density. 
 '''
 pl = point.split() 
 if pl[24] == "Micromegas" : return False 
 '''
 Omgh2obs = float(pl[24])
 if not abs(Omgh2obs - ScanVar.Omgh2) <= ScanVar.Omgh2err : return False 
 '''
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
 Omgh2obs = float(pspl[24])
 
 chi = [(mHiggsobs - ScanVar.mHiggspdg)/(ScanVar.mHiggsSigma), 
         (vEWobs - ScanVar.vEWpdg)/(ScanVar.vEWSigma),
         (Omgh2obs - ScanVar.Omgh2)/(ScanVar.Omgh2Sigma)]
 
 return chi 
 
 
 
def calcChiSqrd(pspoint) : 
 '''
 Reads the parameter space point from a string, 
 calculates the chi squared (Higgs mass and EW VEV), 
 adds it at the end of the string. 
 '''
 chi = calcChiVec(pspoint)
 chisqrd = chi[0]**2 + chi[1]**2 + chi[2]**2
 pspoint += str(format(chisqrd,".8E")) + "\t" 
 
 return pspoint
 
 
 
def runSVM(params, out4lmfit = True) : # 14 parameters 
 '''
 Runs SPheno-Vevacious-Micromegas for a single parameter space point. 
 Returns a string with format:
 "# Lambda1\tLambda2\tLambda3\tLambda4\tLambda5\tLambda6\tLambda7\t" 
     + "Mu122\tTanBeta\tMS2\tLambdaS\tLambda1S\tLambda2S\tLambda12S\t" 
     + "Mhh_1\tMhh_2\tMAh_2\tMHm_2\tMss\t" + "TreeLvlUnit\tTreeLvlUnitwTri\t"
     + "Stability\tFieldValues: v_1\tv_2\tOmega*h^2\tChi squared\n" 
 or the chi vector for lmfit. 
 '''
 makeLesHouchesin(params)
 
 pspoint = runSPheno() 
 if not checkSPout(pspoint) : return np.full(3, 1E+50)
 
 pspoint += runVevacious() 
 if not checkVEVout(pspoint) : return np.full(3, 1E+50) 
 
 pspoint += runMicrOmegas()
 if not checkMicrout(pspoint) : return np.full(3, 1E+50) 
 
 if out4lmfit : return calcChiVec(pspoint) 
 
 pspoint = calcChiSqrd(pspoint) + "\n" 
 
 return [pspoint, calcChiVec(pspoint)] 
 
 
 
def main() : 
 '''
 Nothing. 
 '''
 print("I do nothing...")
 print("Seriously, I do no-thing...")
 return 0  
 
 
 
if __name__ == '__main__' : 
    main()
    
    
    
