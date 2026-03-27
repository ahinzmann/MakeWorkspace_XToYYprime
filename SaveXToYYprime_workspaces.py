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


parser = ArgumentParser()
parser.add_argument('-y', '--years', dest='years', action='store', type=str, choices=['2016', '2016APV', '2017', '2018'], default='2017')
parser.add_argument('-t', '--topology', dest='topology', action='store', type=str, choices=['boosted', 'resolved'], default='resolved')
parser.add_argument('-i', '--inputpath', dest='inputpath', action='store', type=str, default='/data/dust/user/tayu/diboson_output_make3DTemplate_rebinned_MEshower_jun_10p/debug_test/hists/')

args = parser.parse_args()
year = args.years
topology = args.topology
inputpath = args.inputpath

Cut = "Cut 008 HT_cut"
Variable = "Mass_3D"

path_fw = os.environ['CMSSW_BASE']+"/src/MakeWorkspace_XToYYprime/"
store_path = path_fw +"/"+ year + "_" + topology

if __name__ == "__main__":

  os.chdir(path_fw)
  CreatDirectory(store_path)
  os.chdir(store_path)

  openFile = inputpath + Cut + "_" + Variable + ".root" 

  file = ROOT.TFile(openFile)
  for signal_sample in signal_samples:
    signal_path = store_path + "/" + signal_sample
    CreatDirectory(signal_path)
    os.chdir(signal_path)
    Yprime_mass =float( re.search(r'MYprime(\d+)', signal_sample).group(1) )
    print (Yprime_mass)
    Yprime_mass_up   = Yprime_mass*(1+0.3)
    Yprime_mass_down = Yprime_mass*(1-0.3)
    mj2_bins_reduce = filter_list(mj2_bins,Yprime_mass_down,Yprime_mass_up) 
    print (mj2_bins_reduce)
    for sample in samples:
      if "XToYYprime" in sample: 
         if sample != signal_sample:
            continue
      for category in categories: 
        for systematic in systematics:
          dir_list[sample][category][systematic]=[]                    
          #Read the directories
          if sample in ["JetHT"]:
            Read_Hist_Directory(file,sample,systematic,category,dir_list[sample][category][systematic])
          else:
            if sample in ["VV"]: 
              Read_Hist_Directory(file,"ZZ",systematic,category,dir_list[sample][category][systematic])
              Read_Hist_Directory(file,"WZ",systematic,category,dir_list[sample][category][systematic])
              Read_Hist_Directory(file,"WW",systematic,category,dir_list[sample][category][systematic])
            Read_Hist_Directory(file,sample,systematic,category,dir_list[sample][category][systematic])
 
          #Read TT, QCD background and signal hist in direstories
          if dir_list[sample][category][systematic]==[]:
            print (f"{sample} doesn't have category:{category} or systematic:{systematic}, so we skip")
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
 
          ## add mjets and mjets_inverse sys
          if sample == "QCD_madgraph_pythia8" and systematic == "nominal":
             Add_mjets_and_mjetsinverse_sys(sample,category,hist3D_names,hist1D_names)
             Convert_3Dhist_to_1Dhist(sample,category,"mjetsUp",     mj1_bins,mj2_bins_reduce,mjj_bins,hist3D_names,hist_covert3Dto1D)
             Convert_3Dhist_to_1Dhist(sample,category,"mjetsDown",   mj1_bins,mj2_bins_reduce,mjj_bins,hist3D_names,hist_covert3Dto1D)
             Convert_3Dhist_to_1Dhist(sample,category,"mjetsinvUp",  mj1_bins,mj2_bins_reduce,mjj_bins,hist3D_names,hist_covert3Dto1D)
             Convert_3Dhist_to_1Dhist(sample,category,"mjetsinvDown",mj1_bins,mj2_bins_reduce,mjj_bins,hist3D_names,hist_covert3Dto1D)
 
          #convert the 3D hist to 1D
          Convert_3Dhist_to_1Dhist(sample,category,systematic,mj1_bins,mj2_bins_reduce,mjj_bins,hist3D_names,hist_covert3Dto1D)
 
          #Get the projection of 3D histogram
          hist1D_names[sample][category][systematic]["fatjet"] = hist3D_names[sample][category][systematic].ProjectionX("hist_"+sample+"_"+category+"_"+systematic+"_mfatjet")
          hist1D_names[sample][category][systematic]["2jets"] = hist3D_names[sample][category][systematic].ProjectionY("hist_"+sample+"_"+category+"_"+systematic+"_m2jets")
          hist1D_names[sample][category][systematic]["3jets"] = hist3D_names[sample][category][systematic].ProjectionZ("hist_"+sample+"_"+category+"_"+systematic+"_m3jets")
 
 
          if systematic == "nominal" and sample == "JetHT":
            Roodatahist_names[sample][category][systematic] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names[sample][category][systematic])
            Roodatahist_covert3Dto1D[sample][category][systematic] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(x), hist_covert3Dto1D[sample][category][systematic])
            Roodatahist1D_names[sample][category][systematic]["3jets"] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(Mvv), hist1D_names[sample][category][systematic]["3jets"])
          elif systematic == "nominal" and sample != "JetHT" :
            Roodatahist_names[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category,sample+"_"+category,ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names[sample][category][systematic])
            Roodatahist_covert3Dto1D[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category,sample+"_"+category,ROOT.RooArgList(x), hist_covert3Dto1D[sample][category][systematic])
            Roodatahist1D_names[sample][category][systematic]["3jets"] = ROOT.RooDataHist(sample+"_"+category+"_m3jets",sample+"_"+category+"_m3jets",ROOT.RooArgList(Mvv), hist1D_names[sample][category][systematic]["3jets"])
          else:
            Roodatahist_names[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category+"_"+sample+systematic,sample+"_"+category+"_"+sample+systematic,ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names[sample][category][systematic])
            Roodatahist_covert3Dto1D[sample][category][systematic] = ROOT.RooDataHist(sample+"_"+category+"_"+sample+systematic,sample+"_"+category+"_"+sample+systematic,ROOT.RooArgList(x), hist_covert3Dto1D[sample][category][systematic])
            Roodatahist1D_names[sample][category][systematic]["3jets"] = ROOT.RooDataHist(sample+"_"+category+"_m3jets"+"_"+sample+systematic,sample+"_"+category+"_m3jets"+"_"+sample+systematic,ROOT.RooArgList(Mvv), hist1D_names[sample][category][systematic]["3jets"])
           
          # transfer the roodatahist to pdf
          if sample not in ["JetHT"]:
            pdf_names[sample][category][systematic]= ROOT.RooHistPdf(sample+"_"+category+"_"+systematic, sample+"_"+category+"_"+systematic, ROOT.RooArgSet(Mj1, Mj2, Mvv), Roodatahist_names[sample][category][systematic])
            pdf_covert3Dto1D[sample][category][systematic]= ROOT.RooHistPdf(sample+"_"+category+"_"+systematic, sample+"_"+category+"_"+systematic, ROOT.RooArgSet(x), Roodatahist_covert3Dto1D[sample][category][systematic])
 
    import pprint
    #pprint.pprint(result_dict)
    #pprint.pprint(dir_list)
 
    # After getting the different QCD generator and shower hist, then add them in systematics' hist
    Add_generator_shower_sys(categories,hist_covert3Dto1D,hist3D_names,hist1D_names)
 
    #plot the nominal and sys up/down 
    #plot_sys(signal_path,bkg_samples,categories,systematics_names,varname_list,hist1D_names) 
    #plot_sys(signal_path,[signal_sample],categories,systematics_names,varname_list,hist1D_names) 
 
    #make different generator and shower QCD distribution in the same plot
    #plot_QCD_diff_generator_shower(categories,varname_list,hist1D_names) 
 
    #After adding the different QCD generator and shower in systematics' hist, convert them to Roodatahist(including 1D and 3D) 
    for category in categories:
      ###shower
      Roodatahist_names["QCD_madgraph_pythia8"][category]["showerUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_showerUp"  , "QCD_madgraph_pythia8_"+category+"_showerUp"  , ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["showerUp"])
      Roodatahist_names["QCD_madgraph_pythia8"][category]["showerDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_showerDown", "QCD_madgraph_pythia8_"+category+"_showerDown", ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["showerDown"])
      
 
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["showerUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_showerUp"  , "QCD_madgraph_pythia8_"+category+"_showerUp"  , ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["showerUp"])
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["showerDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_showerDown", "QCD_madgraph_pythia8_"+category+"_showerDown", ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["showerDown"])
      ###ME
      Roodatahist_names["QCD_madgraph_pythia8"][category]["MEUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_MEUp"  , "QCD_madgraph_pythia8_"+category+"_MEUp"  , ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["MEUp"])
      Roodatahist_names["QCD_madgraph_pythia8"][category]["MEDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_MEDown", "QCD_madgraph_pythia8_"+category+"_MEDown", ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["MEDown"])
      
 
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_MEUp"  , "QCD_madgraph_pythia8_"+category+"_MEUp"  , ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEUp"])
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_MEDown", "QCD_madgraph_pythia8_"+category+"_MEDown", ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEDown"])
      ###MEshower
      Roodatahist_names["QCD_madgraph_pythia8"][category]["MEshowerUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_MEshowerUp"  , "QCD_madgraph_pythia8_"+category+"_MEshowerUp"  , ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["MEshowerUp"])
      Roodatahist_names["QCD_madgraph_pythia8"][category]["MEshowerDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_MEshowerDown", "QCD_madgraph_pythia8_"+category+"_MEshowerDown", ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["MEshowerDown"])
      
 
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEshowerUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_MEshowerUp"  , "QCD_madgraph_pythia8_"+category+"_MEshowerUp"  , ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEshowerUp"])
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEshowerDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_MEshowerDown", "QCD_madgraph_pythia8_"+category+"_MEshowerDown", ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEshowerDown"])
    #print(hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["herwigUp"])
    
      ###mjets
      Roodatahist_names["QCD_madgraph_pythia8"][category]["mjetsUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_mjetsUp"  , "QCD_madgraph_pythia8_"+category+"_mjetsUp"  , ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["mjetsUp"])
      Roodatahist_names["QCD_madgraph_pythia8"][category]["mjetsDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_mjetsDown", "QCD_madgraph_pythia8_"+category+"_mjetsDown", ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["mjetsDown"])
      
 
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["mjetsUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_mjetsUp"  , "QCD_madgraph_pythia8_"+category+"_mjetsUp"  , ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["mjetsUp"])
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["mjetsDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_mjetsDown", "QCD_madgraph_pythia8_"+category+"_mjetsDown", ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["mjetsDown"])
 
      ###mjets inverse
      Roodatahist_names["QCD_madgraph_pythia8"][category]["mjetsinvUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_mjetsinvUp"  , "QCD_madgraph_pythia8_"+category+"_mjetsinvUp"  , ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["mjetsinvUp"])
      Roodatahist_names["QCD_madgraph_pythia8"][category]["mjetsinvDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_mjetsinvDown", "QCD_madgraph_pythia8_"+category+"_mjetsinvDown", ROOT.RooArgList(Mj1, Mj2, Mvv), hist3D_names["QCD_madgraph_pythia8"][category]["mjetsinvDown"])
      
 
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["mjetsinvUp"]   = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_mjetsinvUp"  , "QCD_madgraph_pythia8_"+category+"_mjetsinvUp"  , ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["mjetsinvUp"])
      Roodatahist_covert3Dto1D["QCD_madgraph_pythia8"][category]["mjetsinvDown"] = ROOT.RooDataHist("QCD_madgraph_pythia8_"+category+"_mjetsinvDown", "QCD_madgraph_pythia8_"+category+"_mjetsinvDown", ROOT.RooArgList(x), hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["mjetsinvDown"])
 
    #Make pesudo data
    bkg_hist_HP_list   = [hist_covert3Dto1D[sample]["HP"]["nominal"] for sample in bkg_samples]
    bkg_hist_LP_list   = [hist_covert3Dto1D[sample]["LP"]["nominal"] for sample in bkg_samples] 
    bkg_hist_rest_list = [hist_covert3Dto1D[sample]["rest"]["nominal"] for sample in bkg_samples] 
 
    hist_pesudo_data_HP   = MakePesudoData_bkgonly(bkg_hist_HP_list)
    hist_pesudo_data_LP   = MakePesudoData_bkgonly(bkg_hist_LP_list)
    hist_pesudo_data_rest = MakePesudoData_bkgonly(bkg_hist_rest_list)
 
 
    #### convert pesudo data hist to Roohist
    Roodatahist_covert3Dto1D["JetHT"]["HP"]["nominal"] = None 
    Roodatahist_covert3Dto1D["JetHT"]["LP"]["nominal"] = None 
    #Roodatahist_covert3Dto1D["JetHT"]["rest"]["nominal"] = None 
 
    Roodatahist_covert3Dto1D["JetHT"]["HP"]["nominal"] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(x), hist_pesudo_data_HP)
    Roodatahist_covert3Dto1D["JetHT"]["LP"]["nominal"] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(x), hist_pesudo_data_LP)
    #Roodatahist_covert3Dto1D["JetHT"]["rest"]["nominal"] = ROOT.RooDataHist("data_obs","data_obs",ROOT.RooArgList(x), hist_pesudo_data_rest)
 
 
    # after HP fit, store the pdf in workspace
    for category in categories: 
      MakeWorkspace(category,signal_sample,bkg_samples,systematics,Roodatahist_covert3Dto1D)  
      WriteDatacard(category,signal_sample,Roodatahist_covert3Dto1D)
 
 
    for category in categories: 
      for sample in bkg_samples:
        #print (sample+" "+category+ " yield is ",  hist3D_names[sample][category]["nominal"].Integral()) 
        #print (sample+" "+category+ " yield is ",  hist1D_names[sample][category]["nominal"]["3jets"].Integral()) 
        #print (sample+" "+category+ " yield is ",  Roodatahist_names[sample][category]["nominal"].sumEntries()) 
        print (sample+" "+category+ " yield is ",  Roodatahist_covert3Dto1D[sample][category]["nominal"].sumEntries()) 
        #print (sample+" "+category+ " yield is ",  hist_covert3Dto1D[sample][category]["nominal"].Integral()) 
 
    print (signal_sample + " and background datacard and workspace have done !")
   
