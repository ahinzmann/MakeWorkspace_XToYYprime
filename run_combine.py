import sys,math,ctypes,array
import ROOT
import os
import re
import matplotlib.pyplot as plt
from ROOT import gROOT, gPad, gStyle
from array import array
#from Save_tools import *
from parameter import signal_samples
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-y', '--years', dest='years', action='store', type=str, choices=['2016', '2016APV', '2017', '2018'], default='2017')
parser.add_argument('-t', '--topology', dest='topology', action='store', type=str, default='all')
parser.add_argument('-rt', '--runtype', dest='runtype', action='store', type=str, choices=['limit', 'fit'], default='limit')

SRcategory_sets=[["2fatjetsHPHP","2fatjetsHPLP"],["3jetsExclHP","3jetsExclLP"],["3jetsExclHP","3jetsExclLP","2fatjetsHPHP","2fatjetsHPLP"]]
CRcategory_sets=[["2fatjetsHPRest"],["3jetsExclRest"],["3jetsExclRest","2fatjetsHPRest"]]

args = parser.parse_args()
year = args.years
topology = args.topology
runtype = args.runtype

path_fw = os.environ['CMSSW_BASE']+"/src/MakeWorkspace_XToYYprime/"
store_path = path_fw +"/"+ year + "_" + topology

def main():
  os.chdir(store_path)
  mass_points = [re.search(r'MX\d+', s).group() for s in signal_samples]
  for i, mass_point in enumerate(mass_points): 
    os.chdir(store_path+"/"+signal_samples[i])
    if runtype == "limit":
     for SRcategories in SRcategory_sets:
      SRname="SR_"+str(SRcategories).strip("[]").replace(", ","_").replace("'","")+f"_{mass_point}"
      os.system("combineCards.py "+str(["c"+cat+"=datacard_"+cat+f"_{mass_point}.txt" for cat in SRcategories]).strip("[]").replace(","," ")+f" > datacard_"+SRname+".txt")
      os.system(f"text2workspace.py datacard_"+SRname+".txt -o workspace_"+SRname+".root")
      os.system(f"combine -M AsymptoticLimits workspace_"+SRname+".root -t -1 > combine_result_"+SRname+".txt")
      os.system(f"cat combine_result_"+SRname+".txt")
    if runtype == "fit":
     for CRcategories in CRcategory_sets:
      CRname="CR_"+str(CRcategories).strip("[]").replace(", ","_").replace("'","")
      print(f"start to run the {signal_samples[i]} FitDiagnostics" )
      os.system("combineCards.py "+str(["c"+cat+"=datacard_"+cat+f"_{mass_point}.txt" for cat in CRcategories]).strip("[]").replace(","," ")+f" > datacard_"+CRname+".txt")
      os.system(f"text2workspace.py datacard_"+CRname+".txt -o workspace_"+CRname+".root")
      os.system(f"combine -M FitDiagnostics workspace_"+CRname+".root --saveShapes --saveWithUncertainties --cminDefaultMinimizerStrategy 0 --ignoreCovWarning -n _"+CRname)
  
if __name__ == "__main__":
  main()
