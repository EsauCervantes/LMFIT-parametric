#!/usr/bin/env python3

# Relic density scan for viable parameter space points 
# scan.dat as argument 
# Run in directory /Vevacious 



import subprocess
import numpy as np 
import scan_spvev as sv 



#------------------Functions---------------------
def makeLesHouchesin(point) : # lmfit.parameter.Parameter object with 14 parameters 
 '''
 Writes new LesHouches.in.2HSDM_low from the original
 with given parameters (list of floats). 
 Sets Lam1S=Lam2S and Lam12S=0 if Lam6=0=Lam7. 
 '''
 pspoint = point.valuesdict() # Dictionary of parameters 
 
 fin = open(sv.ScanPathTo.LHin,"r") 
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
             + " 4   " + str(format(pspoint['Lam4'],".7E")) + "    # Lambda4Input\n" \
             + " 5   " + str(format(pspoint['Lam5'],".7E")) + "    # Lambda5Input\n" \
             + " 6   " + str(format(pspoint['Lam6'],".7E")) + "    # Lambda6Input\n" \
             + " 7   " + str(format(pspoint['Lam7'],".7E")) + "    # Lambda7Input\n" \
             + " 9   " + str(format(pspoint['Mu122'],".7E")) + "    # M122Input\n" \
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
 
 fin = open(sv.ScanPathTo.LHin,"w")
 fin.write(dataout)
 fin.close()
 
 del dataout
 return 0 



def runMicrOmegas() : 
 '''
 Runs MicrOmegas and returns a string which is 
 the resulting relic density read from the file omg.out. 
 '''
 subprocess.run(["cp", sv.ScanPathTo.LHspc, sv.ScanPathTo.MicroDir]) # Copy SPheno output to MicrOmegas directory 
 omgout = ""
 print("Running Micromegas...") 
 
 try : 
     migaspr = subprocess.run([sv.ScanPathTo.Microexe], 
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
     
 fomg = open(sv.ScanPathTo.Omegaout,"r") # Omegaout = "/home/omarset/Vevacious/omg.out"
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



def checkRD(point) : 
 '''
 Checks relic density. 
 '''
 pl = point.split() 
 if pl[24] == "Micromegas" : return False 
 
 Omgh2obs = float(pl[24])
 if not abs(Omgh2obs - sv.ScanVar.Omgh2) <= sv.ScanVar.Omgh2err : return False 
 
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
 
 chi = [(mHiggsobs - sv.ScanVar.mHiggspdg)/(sv.ScanVar.mHiggsSigma), 
         (vEWobs - sv.ScanVar.vEWpdg)/(sv.ScanVar.vEWSigma),
         (Omgh2obs - sv.ScanVar.Omgh2)/(sv.ScanVar.Omgh2Sigma)]
 
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
 
 
 
def runSVM(params) : # 14 parameters 
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
 
 pspoint = sv.runSPheno() 
 if not sv.ViaPointC(pspoint, fullcond = False) : return np.full(3, 1E+50)
 
 pspoint += sv.runVevacious() 
 if not sv.ViaPointC(pspoint, fullcond = True) : return np.full(3, 1E+50) 
 
 pspoint += runMicrOmegas()
 if not checkRD(pspoint) : return np.full(3, 1E+50) 
 #pspoint = calcChiSqrd(pspoint) + "\n" 
 #return [pspoint, calcChiVec(pspoint)] 
 return calcChiVec(pspoint) 
 
 
 
######################## UNDER DEVELOPMENT #################################
'''
Other ideas. 
'''
def makeLHin_TanBeta(TanBetain) :  
 '''
 Writes new LesHouches.in.2HSDM_low changing only tanbeta. 
 '''
 fin = open(sv.ScanPathTo.LHin,"r") 
 dataout = ""
 replaced = False 
 blockMP = False 
 count = 0 

 for line in fin : 
     for index,word in enumerate(line.split()) :         
         if word.startswith("#") : break
         elif word == "MINPAR" :
             blockMP = True 
             break
         if blockMP and index == 0 and word == "10" : 
             dataout += " 10   " + str(format(TanBetain,".7E")) + "    # TanBeta\n" 
             replaced = True
             blockMP = False 
             count += 1
             break
     if not replaced : dataout += line
     else : 
         if count > 1 : 
             dataout += line
             replaced = False 
         count += 1
 fin.close()
    
 fin = open(sv.ScanPathTo.LHin,"w")
 fin.write(dataout)
 fin.close()

 del dataout 
 return 0 
  
 
 
def scanTanBeta_C(TanBetamin, TanBetamax) : # Input consistent with Vevacious minimum 
 '''
 Finds a tanbeta value consistent with Vevacious output. 
 '''
 ####################### Under construction ##############################
 '''
 makeLHin_TanBeta(TanBetamin)
 
 makeLHin_TanBeta(TanBeta)
 svout = sv.runVevacious()
 svl = svout.split() 
 if svl[0] != "stable" : return 0 
 v1 = float(svl[1])
 v2 = float(svl[2])
 TanBetaVout = v2/v1
 dabs = abs(TanBeta - TanBetaVout)
 
 if dabs > sv.ScanVar.tberr : 
     TanBeta += dabs/100
     TanBeta -= dabs/100
 else : return TanBeta 
 
 makeLHin_TanBeta(TanBetamax)
 '''
 return 0 



def makeLHinMicrOmg(MS2in) : # parameter space point, MS2 parameter 
 '''
 Writes new LesHouches.in.2HSDM_low to run MicrOmegas, 
 replaces only MS2 parameter. 
 '''
 fin = open(sv.ScanPathTo.LHin,"r") 
 dataout = ""
 replaced = False 
 blockMP = False 
 count = 0 

 for line in fin : 
     for index,word in enumerate(line.split()) :         
         if word.startswith("#") : break
         elif word == "MINPAR" :
             blockMP = True 
             break
         if blockMP and index == 0 and word == "11" : 
             dataout += " 11   " + str(format(MS2in,".7E")) + "    # MS2Input\n" 
             replaced = True
             blockMP = False 
             count += 1
             break
     if not replaced : dataout += line
     else : 
         if count > 1 : 
             dataout += line
             replaced = False 
         count += 1
 fin.close()
    
 fin = open(sv.ScanPathTo.LHin,"w")
 fin.write(dataout)
 fin.close()

 del dataout 
 
 
 
def scanRelicDens(pspoint, MS2min, MS2max) : 
 '''
 Bisection method on parameter MS2 to obtain 
 the observed dark matter relic density. 
 '''
    
 Nmaxstep = 20 
 a = MS2min
 b = MS2max
 step = 1 
 
 print("Bisection method for DM relic density:")
 
 while step <= Nmaxstep :  
  print("Step " + str(step))
  
  if step == 1 : 
      makeLHinMicrOmg(a)
      if runSPheno() : 
          omga = runMicrOmegas() 
          if abs(sv.ScanVar.Omgh2 - omga ) <= sv.ScanVar.Omgh2err :
              return omga
      else : print("Could not evaluate inferior limit.")
  
  MS2new = (a + b)/2 
  ###################### Under construction #######################
  makeLHinMicrOmg(b)
  if runSPheno() : omgb = runMicrOmegas() 
  else : break # Provisional, need to try it to understand how to deal with errors
  
  makeLHinMicrOmg(MS2new)
  if runSPheno() : omgc = runMicrOmegas()
  else : break # Same as above 
  
  if abs(sv.ScanVar.Omgh2 - omgc) <= sv.ScanVar.Omgh2err : 
      return pointc 
  
  step += 1 
  
  if omgc < sv.ScanVar.Omgh2 < omgb or omgb < sv.ScanVar.Omgh2 < omgc : 
      a = MS2new 
  else : 
      b = MS2new 
  
  if step > Nmaxstep : break 
    
 return 0
 
 
 
def main() : #------------------Main---------------------
 '''
 Implements the scan for each point in a file. 
 '''
 
 '''
 fileobj = open(sys.argv[1], "r")
 
 for index,line in enumerate(fileobj) : 
  if line.startswith("#") : continue 
  print(line.strip("\n")) 
  if index == 1 : break 
 '''
  #makeLesHouchesin(line)
  #scanRelicDens(line,500,10000)
 
 return 0 



if __name__ == '__main__' : 
    main()



