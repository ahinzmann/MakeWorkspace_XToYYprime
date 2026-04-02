import sys,math,ctypes,array
import ROOT
import os
import matplotlib.pyplot as plt
from ROOT import gROOT, gPad, gStyle
from array import array
from Save_tools import *
from parameter import *
from argparse import ArgumentParser

ROOT.gStyle.SetCanvasColor(0)
ROOT.gStyle.SetFrameBorderMode(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTitleX(0.5) #title X location
ROOT.gStyle.SetTitleY(0.96) #title Y location
ROOT.gStyle.SetPaintTextFormat(".2f")
ROOT.gROOT.SetBatch(True)


parser = ArgumentParser()
parser.add_argument('-y', '--years', dest='years', action='store', type=str, choices=['2016', '2016APV', '2017', '2018'], default='2017')
parser.add_argument('-t', '--topology', dest='topology', action='store', type=str, default='all')
parser.add_argument('-i', '--inputpath', dest='inputpath', action='store', type=str, default='/data/dust/user/hinzmann/diboson/pepper/diboson-pepper/diboson_output_make3DTemplate/')

args = parser.parse_args()
year = args.years
topology = args.topology
inputpath = args.inputpath

Cut = "Cut_006_At_least_one_fatjet"
Variables = ["Mass_3D_2fatjets", "Mass_3D_3jets"]

path_fw = os.environ['CMSSW_BASE']+"/src/MakeWorkspace_XToYYprime/"
store_path = path_fw +"/"+ year + "_" + topology

if __name__ == "__main__":

  os.chdir(path_fw)
  CreatDirectory(store_path)
  os.chdir(store_path)

  fs={}
  fs["all"]=["XToYYprime",
           "QCD_madgraph_pythia8",
           "TT",
           "JetHT"]
  fs["all3"]=["VV",
            "ST",
            "ZJets",
            "WJets"]
  fs["all2"]=["QCD_madgraph_herwig7",
            "QCD_herwig7_Pt",
            "QCD_pythia8_Pt"]
  files={}
  for l in fs.keys():
    files[l]=[]
    for v in Variables:
      files[l] += [ROOT.TFile(inputpath + l + "/hists/" + Cut + "_" + v + ".root")]
  for signal_sample in signal_samples:
    signal_path = store_path + "/" + signal_sample
    CreatDirectory(signal_path)
    os.chdir(signal_path)
    Yprime_mass =float( re.search(r'MYprime(\d+)', signal_sample).group(1) )
    print (Yprime_mass)
    Yprime_mass_up   = Yprime_mass*(1+0.3)
    Yprime_mass_down = Yprime_mass*(1-0.3)
    mj2_bins_2fatjets_reduce = filter_list(mj2_bins_2fatjets,Yprime_mass_down,Yprime_mass_up) 
    if len(mj2_bins_2fatjets_reduce)==0:
      mj2_bins_2fatjets_reduce=mj2_bins_2fatjets[:7]
    mj2_bins_3jets_reduce = filter_list(mj2_bins_3jets,Yprime_mass_down,Yprime_mass_up) 
    print (mj2_bins_3jets_reduce,mj2_bins_2fatjets_reduce)
    for sample in samples:
      if "XToYYprime" in sample: 
         if sample != signal_sample:
            continue
      for f in fs.keys():
        if any([(f in sample) for f in fs[f]]):
          filel=f
      for category in categories: 
        if "2fatjets" in category:
          file=files[filel][0]
          mj1_bins=mj1_bins_2fatjets
          mj2_bins_reduce=mj2_bins_2fatjets_reduce
          mjj_bins=mjj_bins_2fatjets
          x=x_2fatjets
        else:
          file=files[filel][1]
          mj1_bins=mj1_bins_3jets
          mj2_bins_reduce=mj2_bins_3jets_reduce
          mjj_bins=mjj_bins_3jets
          x=x_3jets
        for systematic in systematics:
          dir_list[sample][category][systematic]=[]                    
          #Read the directories
          if "VV" in sample: 
              Read_Hist_Directory(file,"ZZ",systematic.replace("Up","_up").replace("Down","_down"),category,dir_list[sample][category][systematic])
              Read_Hist_Directory(file,"WZ",systematic.replace("Up","_up").replace("Down","_down"),category,dir_list[sample][category][systematic])
              Read_Hist_Directory(file,"WW",systematic.replace("Up","_up").replace("Down","_down"),category,dir_list[sample][category][systematic])
          else:
              Read_Hist_Directory(file,sample,systematic.replace("Up","_up").replace("Down","_down"),category,dir_list[sample][category][systematic])
 
          #Read TT, QCD background and signal hist in direstories
          if dir_list[sample][category][systematic]==[]:
            print (f"{sample} doesn't have category:{category} or systematic:{systematic}, so we skip")
            print(sample,filel,category)
            continue
          else:
            print (f"Now running the {sample}_{category}_{systematic}")
            hist3D_names[sample][category][systematic] = Read_3DHist(dir_list[sample][category][systematic],sample)
          #if sample == "XToYYprime_MX3000":
          #  hist3D_names[sample][category][systematic].Scale(0.000415) 
          if sample == "QCD_madgraph_herwig7":
            hist3D_names[sample][category][systematic].Scale(1500)
          if sample == "QCD_herwig7_Pt":
            hist3D_names[sample][category][systematic].Scale(1.12)
          if sample == "QCD_pythia8_Pt":
            hist3D_names[sample][category][systematic].Scale(1)
 
          ## add functional sys
          if sample == "QCD_madgraph_pythia8" and systematic == "nominal":
             Add_functional_sys(sample,category,hist3D_names,hist1D_names,varname_list)
             for vn in varname_list:
              for dn in ["Up","Down","invUp","invDown"]:
               Convert_3Dhist_to_1Dhist(sample,category,"m"+vn+dn,mj1_bins,mj2_bins_reduce,mjj_bins,hist3D_names,hist_covert3Dto1D)
 
          #convert the 3D hist to 1D
          Convert_3Dhist_to_1Dhist(sample,category,systematic,mj1_bins,mj2_bins_reduce,mjj_bins,hist3D_names,hist_covert3Dto1D)
 
         #Get the projection of 3D histogram
          hist1D_names[sample][category][systematic]["fatjet"] = hist3D_names[sample][category][systematic].ProjectionX("hist_"+sample+"_"+category+"_"+systematic+"_mfatjet")
          hist1D_names[sample][category][systematic]["2jets"] = hist3D_names[sample][category][systematic].ProjectionY("hist_"+sample+"_"+category+"_"+systematic+"_m2jets")
          hist1D_names[sample][category][systematic]["3jets"] = hist3D_names[sample][category][systematic].ProjectionZ("hist_"+sample+"_"+category+"_"+systematic+"_m3jets")
 
          if systematic == "nominal" and sample == "JetHT": # data
            #Roodatahist_names[sample][category][systematic] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names[sample][category][systematic])
            #Roodatahist_covert3Dto1D[sample][category][systematic] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(x), hist_covert3Dto1D[sample][category][systematic])
            hist_covert3Dto1D[sample][category][systematic].SetName("data_obs_"+category)
            #Roodatahist1D_names[sample][category][systematic]["3jets"] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(Mvv), hist1D_names[sample][category][systematic]["3jets"])
          elif systematic == "nominal" and sample != "JetHT": # nominal
            #Roodatahist_names[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category,sample+"_"+category,ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names[sample][category][systematic])
            #Roodatahist_covert3Dto1D[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category,sample+"_"+category,ROOT.RooArgList(x), hist_covert3Dto1D[sample][category][systematic])
            hist_covert3Dto1D[sample][category][systematic].SetName(sample+"_"+category)
            #Roodatahist1D_names[sample][category][systematic]["3jets"] = ROOT.RooDataHist(sample+"_"+category+"_m3jets",sample+"_"+category+"_m3jets",ROOT.RooArgList(Mvv), hist1D_names[sample][category][systematic]["3jets"])
          elif any([(f in systematic) for f in ["JuncTotal","PSfsr"]]): # correlated among processes
            #Roodatahist_names[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category+"_"+systematic,sample+"_"+category+"_"+systematic,ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names[sample][category][systematic])
            #Roodatahist_covert3Dto1D[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category+"_"+systematic,sample+"_"+category+"_"+systematic,ROOT.RooArgList(x), hist_covert3Dto1D[sample][category][systematic])
            hist_covert3Dto1D[sample][category][systematic].SetName(sample+"_"+category+"_"+systematic)
            #Roodatahist1D_names[sample][category][systematic]["3jets"] = ROOT.RooDataHist(sample+"_"+category+"_m3jets"+"_"+systematic,sample+"_"+category+"_m3jets"+"_"+systematic,ROOT.RooArgList(Mvv), hist1D_names[sample][category][systematic]["3jets"])
          else: # uncorrelated among processes
            #Roodatahist_names[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category+"_"+sample+systematic,sample+"_"+category+"_"+sample+systematic,ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names[sample][category][systematic])
            #Roodatahist_covert3Dto1D[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category+"_"+sample+systematic,sample+"_"+category+"_"+sample+systematic,ROOT.RooArgList(x), hist_covert3Dto1D[sample][category][systematic])
            hist_covert3Dto1D[sample][category][systematic].SetName(sample+"_"+category+"_"+sample+systematic)
            #Roodatahist1D_names[sample][category][systematic]["3jets"] = ROOT.RooDataHist(sample+"_"+category+"_m3jets"+"_"+sample+systematic,sample+"_"+category+"_m3jets"+"_"+sample+systematic,ROOT.RooArgList(Mvv), hist1D_names[sample][category][systematic]["3jets"])

          # transfer the roodatahist to pdf (only for 3D limits)
          #if sample not in ["JetHT"]:
          #  pdf_names[sample][category][systematic]= ROOT.RooHistPdf(sample+"_"+category+"_"+systematic, sample+"_"+category+"_"+systematic, ROOT.RooArgSet(Mj1, Mj2, Mvv), Roodatahist_names[sample][category][systematic])
          #  pdf_covert3Dto1D[sample][category][systematic]= ROOT.RooHistPdf(sample+"_"+category+"_"+systematic, sample+"_"+category+"_"+systematic, ROOT.RooArgSet(x), Roodatahist_covert3Dto1D[sample][category][systematic])
          #print("store",category,mj1_bins,mj2_bins_reduce,mjj_bins,hist_covert3Dto1D[sample][category][systematic].GetNbinsX())
          #print("store",category,Roodatahist_covert3Dto1D[sample][category][systematic].numEntries())

    import pprint
    #pprint.pprint(result_dict)
    #pprint.pprint(dir_list)
 
    # After getting the different QCD generator and shower hist, then add them in systematics' hist
    Add_generator_shower_sys(categories,hist_covert3Dto1D,hist3D_names,hist1D_names)
 
    #plot the nominal and sys up/down 
    plot_sys(signal_path,bkg_samples,categories,systematics_names,varname_list,hist1D_names) 
    plot_sys(signal_path,[signal_sample],categories,systematics_names,varname_list,hist1D_names) 
 
    #make different generator and shower QCD distribution in the same plot
    plot_QCD_diff_generator_shower(categories,varname_list,hist1D_names) 
 
    #After adding the different QCD generator and shower in systematics' hist, convert them to Roodatahist(including 1D and 3D) 
    for category in categories:
      if "2fatjets" in category:
          x=x_2fatjets
      else:
          x=x_3jets
      
      for sys in ["shower","ME","MEshower"]+["m"+vn for vn in varname_list]+["m"+vn+"inv" for vn in varname_list]:
        for dn in ["Up","Down"]:
          #Roodatahist_names["QCD_madgraph_pythia8"][category][sys+dn]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_"+sys+dn  , "QCD_madgraph_pythia8_"+category+"_"+sys+dn  , ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category][sys+dn])
          #Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category][sys+dn]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_"+sys+dn  , "QCD_madgraph_pythia8_"+category+"_"+sys+dn  , ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category][sys+dn])
          hist_covert3Dto1D["QCD_madgraph_pythia8"][category][sys+dn].SetName("QCD_madgraph_pythia8_"+category+"_"+sys+dn)
 
    # Make pesudo data
    for category in categories:
      bkg_hist_list   = [hist_covert3Dto1D[sample][category]["nominal"] for sample in bkg_samples]
      hist_pesudo_data   = MakePesudoData_bkgonly(bkg_hist_list)
 
      #### convert pesudo data hist to Roohist
      if not "Rest" in category:
        #Roodatahist_covert3Dto1D["JetHT"][category]["nominal"] = None 
        #Roodatahist_covert3Dto1D["JetHT"][category]["nominal"] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(x), hist_pesudo_data)
        hist_covert3Dto1D["JetHT"][category]["nominal"] = hist_pesudo_data
        hist_covert3Dto1D["JetHT"][category]["nominal"].SetName("data_obs_"+category)

      # after HP fit, store the pdf in workspace
      #MakeWorkspace(category,signal_sample,bkg_samples,systematics,Roodatahist_covert3Dto1D)  
      MakeWorkspaceSimple(category,signal_sample,bkg_samples,systematics,hist_covert3Dto1D)
      WriteDatacard(category,signal_sample,hist_covert3Dto1D)
      #print("store",category,Roodatahist_covert3Dto1D["JetHT"][category]["nominal"].numEntries())
      print("store",category,hist_covert3Dto1D["JetHT"][category]["nominal"].Integral())
 
      for sample in bkg_samples:
        #print (sample+" "+category+ " yield is ",  hist3D_names[sample][category]["nominal"].Integral()) 
        #print (sample+" "+category+ " yield is ",  hist1D_names[sample][category]["nominal"]["3jets"].Integral()) 
        #print (sample+" "+category+ " yield is ",  Roodatahist_names[sample][category]["nominal"].sumEntries()) 
        #print (sample+" "+category+ " yield is ",  Roodatahist_covert3Dto1D[sample][category]["nominal"].sumEntries()) 
        print (sample+" "+category+ " yield is ",  hist_covert3Dto1D[sample][category]["nominal"].Integral())
        #print (sample+" "+category+ " yield is ",  hist_covert3Dto1D[sample][category]["nominal"].Integral()) 
 
    print (signal_sample + " and background datacard and workspace have done !")

