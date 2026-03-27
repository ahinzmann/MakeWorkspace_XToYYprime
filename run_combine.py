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
parser.add_argument('-t', '--topology', dest='topology', action='store', type=str, choices=['boosted', 'resolved'], default='resolved')
parser.add_argument('-rt', '--runtype', dest='runtype', action='store', type=str, choices=['get_limit', 'estimate_bkg'], default='estimate_bkg')


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
    if runtype == "get_limit":
      os.system(f"combineCards.py HP=datacard_HP_{mass_point}.txt LP=datacard_LP_{mass_point}.txt > datacard_SR_{mass_point}.txt")
      os.system(f"text2workspace.py datacard_SR_{mass_point}.txt -o workspace_SR_{mass_point}.root")
      os.system(f"combine -M AsymptoticLimits workspace_SR_{mass_point}.root -t -1 > combine_result_{mass_point}.txt")
      os.system(f"cat combine_result_{mass_point}.txt")
    if runtype == "estimate_bkg":
      print(f"start to run the {signal_samples[i]} FitDiagnostics" )
      os.system(f"text2workspace.py datacard_rest_{mass_point}.txt -o workspace_CR.root")
      os.system("combine -M FitDiagnostics workspace_CR.root --saveShapes --saveWithUncertainties --cminDefaultMinimizerStrategy 0 --ignoreCovWarning -n postfit ")
  
if __name__ == "__main__":
  main()
