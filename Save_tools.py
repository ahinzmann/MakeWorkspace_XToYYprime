import sys,math,ctypes,array
import ROOT
import os
import re
import matplotlib.pyplot as plt
from ROOT import gROOT, gPad, gStyle
from ROOT import TH1F, TH3F, TCanvas
from array import array
from ROOT import RooRealVar, RooDataHist, RooArgList, RooArgSet, RooGaussian, RooHistPdf, RooFit, RooAddPdf

ROOT.gStyle.SetCanvasColor(0)
ROOT.gStyle.SetFrameBorderMode(0)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPalette(1,0,1.)
ROOT.gStyle.SetTitleX(0.5) #title X location
ROOT.gStyle.SetTitleY(0.96) #title Y location
ROOT.gStyle.SetPaintTextFormat(".2f")



def CreatDirectory(path):
    if not os.path.exists(path):
     os.makedirs(path)
     print("%s created successfully" % path)
    else:
     print("%s already exist" % path)


def filter_list(sorted_list, lower_bound, upper_bound):
    """
    select the element from the ordered list in [lower_bound, upper_bound] 

    parameter:
        sorted_list: ordered list
        lower_bound: up limit(include)
        upper_bound: down limit(include)
        
    return:
        selected list
    """
    filtered = [x for x in sorted_list if lower_bound <= x <= upper_bound]
    return filtered


def Read_Hist_Directory(file,samplename,systematicname,categoryname,dir_list):   
   for key_sample in file.GetListOfKeys():
    if samplename in key_sample.GetName():
      obj_sample = key_sample.ReadObj()
      if isinstance(obj_sample, ROOT.TDirectoryFile):
        for key_systematic in obj_sample.GetListOfKeys():
          if systematicname in key_systematic.GetName():
            obj_systematic = key_systematic.ReadObj() 
            if isinstance(obj_systematic, ROOT.TDirectoryFile):
              for key_category in obj_systematic.GetListOfKeys():
                if categoryname in key_category.GetName():
                  obj_category = key_category.ReadObj()
                  if isinstance(obj_category, ROOT.TDirectoryFile):
                    dir_list.append(obj_category)

def Read_Hist_Directory_nosys(file,samplename,systematicname,categoryname,dir_list):   
   for key_sample in file.GetListOfKeys():
    if samplename in key_sample.GetName():
      obj_sample = key_sample.ReadObj()
      if isinstance(obj_sample, ROOT.TDirectoryFile):
        #for key_systematic in obj_sample.GetListOfKeys():
          #if systematicname in key_systematic.GetName():
            #obj_systematic = key_systematic.ReadObj() 
            #if isinstance(obj_systematic, ROOT.TDirectoryFile):
              for key_category in obj_sample.GetListOfKeys():
                if categoryname in key_category.GetName():
                  obj_category = key_category.ReadObj()
                  if isinstance(obj_category, ROOT.TDirectoryFile):
                    dir_list.append(obj_category)

def Read_3DHist(dir_list,samplename):
   for i, directory in enumerate(dir_list):
    #print (dir_list)
    directory.cd()
    hist_list = []
    for key in directory.GetListOfKeys():
      obj = key.ReadObj()
      if isinstance(obj, ROOT.TH3):
        hist_list.append(obj)

    for j, hist in enumerate(hist_list):
      if samplename == "QCD":
        hist.SetLineColor(3)
        hist.SetFillColor(3)
      if samplename == "TT":
        hist.SetLineColor(4)
        hist.SetFillColor(4)
      if samplename == "signal":
        hist.SetLineColor(2)
      if j == 0 and i == 0:
        hist_total = hist.Clone("hist_total")
      else:
        hist_total.Add(hist)
   return hist_total  

def Read_1DHist(dir_list,samplename):
   for i, directory in enumerate(dir_list):
    print (dir_list)
    directory.cd()
    hist_list = []
    for key in directory.GetListOfKeys():
      obj = key.ReadObj()
      if isinstance(obj, ROOT.TH1):
        hist_list.append(obj)

    for j, hist in enumerate(hist_list):
      if j == 0 and i == 0:
        hist_total = hist.Clone("hist_total")
      else:
        hist_total.Add(hist)
   return hist_total  

def smooth_th3_steep(h3):
    """
    Smoothes a TH3 histogram in log-space to preserve 
    steeply falling spectral shapes.
    """
    import numpy as np
    from scipy.ndimage import gaussian_filter
    nx, ny, nz = h3.GetNbinsX(), h3.GetNbinsY(), h3.GetNbinsZ()
    
    # 1. Extract data to a NumPy array (shape: nx+2, ny+2, nz+2 includes under/overflow)
    data = np.zeros((nx + 2, ny + 2, nz + 2))
    for i in range(nx + 2):
        for j in range(ny + 2):
            for k in range(nz + 2):
                data[i, j, k] = h3.GetBinContent(i, j, k)

    # 2. Log-transform (add small epsilon to avoid log(0))
    # This ensures we smooth the "slope" rather than the absolute counts
    #epsilon = 1e-10
    #log_data = np.log(data + epsilon)

    # 3. Apply Gaussian Filter (sigma controls smoothing strength)
    # sigma=0.5 to 1.0 is usually enough for minor denoising
    #smoothed_log = gaussian_filter(log_data, sigma=0.8)

    # 4. Transform back to linear space
    #smoothed_data = np.exp(smoothed_log)

    smoothed_data = gaussian_filter(data, sigma=0.8)

    # 5. Create new TH3 and fill with smoothed values
    h_smooth = h3.Clone(h3.GetName() + "_smooth")
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            for k in range(1, nz + 1):
                h_smooth.SetBinContent(i, j, k, smoothed_data[i, j, k])
    
    return h_smooth

def Convert_3Dhist_to_1Dhist(sample, category, systematic, mj1_bins, mj2_bins, mjj_bins, hist3D_names, hist_covert3Dto1D):
  h3 = hist3D_names[sample][category][systematic]
  #if not "JetHT" in sample: # Smoothen MC
  #  h3=smooth_th3_steep(h3)
  projX = h3.ProjectionX("projX")
  projY = h3.ProjectionY("projY") 
  projZ = h3.ProjectionZ("projZ")

  min_bin_X    = projX.FindBin(mj1_bins[0]) #min value
  max_bin_X    = projX.FindBin(mj1_bins[-1]-1) #max value-1
  n_bins_X     = max_bin_X-min_bin_X+1
  low_range_X  = projX.GetXaxis().GetBinLowEdge(min_bin_X)
  high_range_X = projX.GetXaxis().GetBinUpEdge(max_bin_X)

  min_bin_Y    = projY.FindBin(mj2_bins[0] )
  max_bin_Y    = projY.FindBin(mj2_bins[-1]-1)
  n_bins_Y     = max_bin_Y-min_bin_Y+1
  low_range_Y  = projY.GetXaxis().GetBinLowEdge(min_bin_Y)
  high_range_Y = projY.GetXaxis().GetBinUpEdge(max_bin_Y)

  min_bin_Z    = projZ.FindBin(mjj_bins[0])
  max_bin_Z    = projZ.FindBin(mjj_bins[-1]-1)
  n_bins_Z     = max_bin_Z-min_bin_Z+1
  low_range_Z  = projZ.GetXaxis().GetBinLowEdge(min_bin_Z)
  high_range_Z = projZ.GetXaxis().GetBinUpEdge(max_bin_Z)
 
  #print(f"mj1 bin is from {min_bin_X} to {max_bin_X}, range is from {low_range_X} to {high_range_X} ,bin number is {n_bins_X}")
  #print(f"mj2 bin is from {min_bin_Y} to {max_bin_Y}, range is from {low_range_Y} to {high_range_Y} ,bin number is {n_bins_Y}")
  #print(f"mjj bin is from {min_bin_Z} to {max_bin_Z}, range is from {low_range_Z} to {high_range_Z} ,bin number is {n_bins_Z}")

  #nBinsX = projX.GetNbinsX()
  #nBinsY = projY.GetNbinsY()  
  #nBinsZ = projZ.GetNbinsZ()
  
  #totalBins = nBinsX * nBinsY * nBinsZ
  totalBins = n_bins_X * n_bins_Y * n_bins_Z
 

  #print("the number of total bin is ", n_bins_X , n_bins_Y , n_bins_Z, totalBins)

  longHist = TH1F("longHist_"+sample+"_"+category+"_"+systematic, "Combined Projections;Bin Index;Entries", totalBins, 0, totalBins)
 
  for i in range(min_bin_X, max_bin_X + 1):
    for j in range(min_bin_Y, max_bin_Y + 1):
      for k in range(min_bin_Z, max_bin_Z + 1):
        index=k-min_bin_Z+1 + n_bins_Z*( (j-min_bin_Y+1 -1) + n_bins_Y * (i-min_bin_X+1 -1) )
        longHist.SetBinContent(index, h3.GetBinContent(i,j,k))
        #print(i,j,k,index,h3.GetBinContent(i,j,k))

  '''
  
  CreatDirectory("./bin_content")
  with open(f"bin_content/{sample}_{category}_{systematic}.txt","w") as f:
    for i in range(min_bin_X, max_bin_X + 1):
      for j in range(min_bin_Y, max_bin_Y + 1):
        for k in range(min_bin_Z, max_bin_Z + 1):
          longHist.SetBinContent(k-min_bin_Z+1 + n_bins_Z*( (j-min_bin_Y+1 -1) + n_bins_Y * (i-min_bin_X+1 -1) ), h3.GetBinContent(i,j,k))
          if sample == "QCD_madgraph_pythia8":
 
            print (f"1D long hist bin {k-min_bin_Z+1 + n_bins_Z*( (j-min_bin_Y+1 -1) + n_bins_Y * (i-min_bin_X+1 -1) )} content is {h3.GetBinContent(i,j,k)}\n",file=f )
  '''

  hist_covert3Dto1D[sample][category][systematic] = longHist

  #c1 = ROOT.TCanvas("c1", "Combined Projections", 1200, 600)
  #longHist.Draw()

  #c1.Draw()
  #c1.SaveAs("convert_3Dto1Dplots/"+sample+"_"+category+"_"+systematic+".pdf")


def RooFit_1D(cut,samplename,category,systematic,fitRange,var,varname,pdf,hist):
   roohist = ROOT.RooDataHist(samplename+"_"+category+"_"+systematic+"_"+varname+"_datahist",samplename+"_"+category+"_"+varname+"_datahist", ROOT.RooArgList(var),hist)
   canvas = ROOT.TCanvas()
   canvas.SetLogy()
   frame = var.frame()
   #pdf.fitTo(roohist, ROOT.SumW2Error(kTRUE))
   pdf.fitTo(roohist,ROOT.RooFit.Range(fitRange))
   roohist.plotOn(frame)
   #pdf.plotOn(frame)
   pdf.plotOn(frame,ROOT.RooFit.LineColor(hist.GetLineColor()))
   #pdf.plotOn(frame,ROOT.RooFit.Components("nonres_erfExp"),ROOT.RooFit.LineColor(2))
   temp = hist.Clone("temp")
   temp.SetLineColor(1)
   temp.SetMarkerStyle(20)
   temp.SetMarkerSize(1.3)

   # Get the maxmum of y axis
   max_bin = hist.GetMaximumBin()
   max_value = hist.GetBinContent(max_bin)

   frame.SetTitle("")
   frame.GetXaxis().SetTitle(varname +" mass [GeV]")
   frame.GetXaxis().SetTitleOffset(0.9)
   frame.GetXaxis().SetTitleSize(0.050)
   frame.GetXaxis().SetLabelSize(0.050)
   frame.GetYaxis().SetTitle("Events")
   frame.GetYaxis().SetLabelSize(0.050)
   frame.GetYaxis().SetTitleSize(0.050)
   frame.GetYaxis().SetTitleOffset(0.9)
   frame.SetMaximum(max_value*10)
   frame.SetMinimum(0.0001)
   frame.Draw()
   gPad.RedrawAxis()

   legend = ROOT.TLegend(0.65, 0.75, 0.88, 0.89)
   legend.SetTextSize(0.050)
   legend.SetTextFont(62)
   legend.SetFillColor(0)
   legend.AddEntry(temp,samplename, "PE")
   legend.AddEntry(hist,samplename+" fit", "l")
   legend.SetBorderSize(0)
   legend.Draw()

   canvas.Draw()

   canvas.SaveAs("RooFit1D_plots/"+samplename+"_"+category+"_"+systematic+"_"+varname+"_mass_projection.pdf")

def RooFit_3D(cut,samplename,category,fitRange,var_list,varname_list,pdf,hist):
   roohist = ROOT.RooDataHist(samplename+"_"+category+"_3Ddatahist",samplename+"_"+category+"_3Ddatahist", ROOT.RooArgList(var_list[0],var_list[1],var_list[2]),hist)
   pdf.fitTo(roohist,ROOT.RooFit.Range(fitRange)) 
   # plot the projection of 3 dimension
   for i in range(len(var_list)): 
     plot_pdf_and_projection(cut,samplename,category,var_list[i],varname_list[i],pdf,hist,roohist) 
     #plot_projection(cut,samplename,category,var_list[i],varname_list[i],hist,roohist) 
   del roohist

def Post_Fit_b(Cut,category,var_list,bkg_samples,samples_color,hist1D_names,pdf_names,hist_pesudo_data):
   #make model b
   n_QCD   = ROOT.RooRealVar("n_QCD"  ,"n_QCD"  , 3924694, 0, 10e7)#, 0, 10e8
   n_ttbar = ROOT.RooRealVar("n_ttbar","n_ttbar", 28648   )#, 0, 20e5
   n_st    = ROOT.RooRealVar("n_st"   ,"n_st"   , 4416  )#, 0, 10e5
   n_zjet  = ROOT.RooRealVar("n_zjet" ,"n_zjet" , 31957   )#, 0, 20e5
   n_wjet  = ROOT.RooRealVar("n_wjet" ,"n_wjet" , 60824   )#, 0, 30e5
   n_vv    = ROOT.RooRealVar("n_vv"   ,"n_vv"   , 874     )#, 0, 10e4
   model_b = ROOT.RooAddPdf("model_b","model_b",ROOT.RooArgList(pdf_names["QCD"][category]["nominal"],pdf_names["TT"][category]["nominal"],pdf_names["ST"][category]["nominal"],pdf_names["ZJets"][category]["nominal"],pdf_names["WJets"][category]["nominal"],pdf_names["VV"][category]["nominal"]),ROOT.RooArgList(n_QCD,n_ttbar,n_st,n_zjet,n_wjet,n_vv))

   #print (PDF_list[1].GetName())
   ##The pesudo data
   hist_pesudo_data.Scale(1) 
   roohist_pesudodata = ROOT.RooDataHist(category+"_3Ddatahist",category+"_3Ddatahist", ROOT.RooArgList(var_list[0],var_list[1],var_list[2]),hist_pesudo_data)

   print(roohist_pesudodata.sumEntries())
   ## signal
   #hist_signal_clone = hist_signal.Clone("hist_signal_clone")
   #hist_signal_clone.Scale(50) 
   #roohist_signal = ROOT.RooDataHist(category+"_3Dsignal",category+"_3Dsignal", ROOT.RooArgList(var_list[0],var_list[1],var_list[2]),hist_signal_clone)

   #Get hist min and max value
   min_mj1 = var_list[0].getMin()   
   max_mj1 = var_list[0].getMax()
   min_mj2 = var_list[1].getMin()   
   max_mj2 = var_list[1].getMax()
   min_mjj = var_list[2].getMin()   
   max_mjj = var_list[2].getMax()

   for sample in bkg_samples:
     hist1D_names[sample][category]["nominal"]["mjj"].SetLineColor(samples_color[sample])
   hist1D_names["JetHT"][category]["nominal"]["mjj"].SetLineColor(1)
   hist1D_names["JetHT"][category]["nominal"]["mjj"].SetMarkerStyle(20)
   hist1D_names["JetHT"][category]["nominal"]["mjj"].SetMarkerSize(1.3)

   #creat the projection
   mj1_frame = var_list[0].frame()
   mj2_frame = var_list[1].frame()
   mjj_frame = var_list[2].frame()

   canvas = ROOT.TCanvas("canvas", "Stacked PDF Projection", 800, 600)
   canvas.SetLogy()
   #fit the pesudodata
   model_b.fitTo(roohist_pesudodata)

   canvas.Divide(1,2);
   for i in range(3):
     if i == 0:
       ##plot the data
       roohist_pesudodata.plotOn(mj1_frame)    
       ## plot the total pdf
       model_b.plotOn(mj1_frame,ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.LineWidth(1))
       ##plot the divid hist
       hpull = mj1_frame.pullHist() 
       mj1_Pull_Frame = var_list[i].frame()
       mj1_Pull_Frame.addPlotable(hpull,"P")

       ## plot the components pdf
       for sample in bkg_samples:
         model_b.plotOn(mj1_frame,ROOT.RooFit.Components(pdf_names[sample][category]["nominal"].GetName()),ROOT.RooFit.LineColor(samples_color[sample]),ROOT.RooFit.LineWidth(1))
       #roohist_signal.plotOn(mj1_frame, ROOT.RooFit.DrawOption("P"), ROOT.RooFit.LineWidth(1),ROOT.RooFit.LineColor(13),ROOT.RooFit.MarkerColor(13), ROOT.RooFit.LineStyle(2))

       canvas.cd(1).SetPad(0.01, 0.23, 0.99, 0.98)
       gPad.SetLeftMargin(0.1) 
       gPad.SetLogy()

       # Get the maxmum of y axis
       max_bin   = hist1D_names["JetHT"][category]["nominal"]["mj1"].GetMaximumBin()
       max_value = hist1D_names["JetHT"][category]["nominal"]["mj1"].GetBinContent(max_bin)
       print (max_value)

       # Draw the frame
       mj1_frame.SetTitle("")
       mj1_frame.GetXaxis().SetTitleOffset(0.9)
       mj1_frame.GetXaxis().SetTitleSize(0.050)
       mj1_frame.GetXaxis().SetLabelSize(0)
       #frame.GetYaxis().SetTitle("Events")
       mj1_frame.GetYaxis().SetLabelSize(0.050)
       mj1_frame.GetYaxis().SetTitleSize(0.050)
       mj1_frame.GetYaxis().SetTitleOffset(0.9)
       mj1_frame.SetMaximum(max_value*5)
       mj1_frame.SetMinimum(1)
       mj1_frame.Draw()
   
       canvas.cd(2).SetPad(0.01, 0.03, 0.99, 0.29)
       gPad.SetLeftMargin(0.1)
       gPad.SetBottomMargin(0.27)
   
       mj1_Pull_Frame.SetTitle("")
       mj1_Pull_Frame.GetXaxis().SetTitle("Mj1 [GeV]")
       mj1_Pull_Frame.GetXaxis().SetTitleOffset(0.9)
       mj1_Pull_Frame.GetXaxis().SetTitleSize(0.15)
       mj1_Pull_Frame.GetXaxis().SetLabelSize(0.15)
       #frame.GetYaxis().SetTitle("Events")
       mj1_Pull_Frame.GetYaxis().SetTitle("Pull")
       mj1_Pull_Frame.GetYaxis().SetLabelSize(0.150)
       mj1_Pull_Frame.GetYaxis().SetTitleSize(0.150)
       mj1_Pull_Frame.GetYaxis().SetTitleOffset(0.9)
       mj1_Pull_Frame.Draw()

       line = ROOT.TLine(min_mj1, 0, max_mj1, 0)
       line.SetLineColor(2)
       line.Draw("same")
     if i == 1:
       ##plot the data
       roohist_pesudodata.plotOn(mj2_frame)    
       ## plot the total pdf
       model_b.plotOn(mj2_frame,ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.LineWidth(1))
       ##plot the divid hist
       hpull = mj2_frame.pullHist() 
       mj2_Pull_Frame = var_list[i].frame()
       mj2_Pull_Frame.addPlotable(hpull,"P")

       ## plot the components pdf
       for sample in bkg_samples:
         model_b.plotOn(mj2_frame,ROOT.RooFit.Components(pdf_names[sample][category]["nominal"].GetName()),ROOT.RooFit.LineColor(samples_color[sample]),ROOT.RooFit.LineWidth(1))

       #roohist_signal.plotOn(mj2_frame, ROOT.RooFit.DrawOption("P"), ROOT.RooFit.LineWidth(1),ROOT.RooFit.LineColor(13),ROOT.RooFit.MarkerColor(13), ROOT.RooFit.LineStyle(2))

       canvas.cd(1).SetPad(0.01, 0.23, 0.99, 0.98)
       gPad.SetLeftMargin(0.1) 
       gPad.SetLogy()

       # Get the maxmum of y axis
       max_bin   = hist1D_names["JetHT"][category]["nominal"]["mj2"].GetMaximumBin()
       max_value = hist1D_names["JetHT"][category]["nominal"]["mj2"].GetBinContent(max_bin)
       print (max_value)

       # Draw the frame
       mj2_frame.SetTitle("")
       mj2_frame.GetXaxis().SetTitleOffset(0.9)
       mj2_frame.GetXaxis().SetTitleSize(0.050)
       mj2_frame.GetXaxis().SetLabelSize(0)
       #frame.GetYaxis().SetTitle("Events")
       mj2_frame.GetYaxis().SetLabelSize(0.050)
       mj2_frame.GetYaxis().SetTitleSize(0.050)
       mj2_frame.GetYaxis().SetTitleOffset(0.9)
       mj2_frame.SetMaximum(max_value*10)
       mj2_frame.SetMinimum(1)
       mj2_frame.Draw()
   
       canvas.cd(2).SetPad(0.01, 0.03, 0.99, 0.29)
       gPad.SetLeftMargin(0.1)
       gPad.SetBottomMargin(0.27)
   
       mj2_Pull_Frame.SetTitle("")
       mj2_Pull_Frame.GetXaxis().SetTitle("Mj2 [GeV]")
       mj2_Pull_Frame.GetXaxis().SetTitleOffset(0.9)
       mj2_Pull_Frame.GetXaxis().SetTitleSize(0.15)
       mj2_Pull_Frame.GetXaxis().SetLabelSize(0.15)
       #frame.GetYaxis().SetTitle("Events")
       mj2_Pull_Frame.GetYaxis().SetTitle("Pull")
       mj2_Pull_Frame.GetYaxis().SetLabelSize(0.150)
       mj2_Pull_Frame.GetYaxis().SetTitleSize(0.150)
       mj2_Pull_Frame.GetYaxis().SetTitleOffset(0.9)
       mj2_Pull_Frame.Draw()

       line = ROOT.TLine(min_mj2, 0, max_mj2, 0)
       line.SetLineColor(2)
       line.Draw("same")
     if i == 2:
       var_list[0].setRange("signal", 50, 500)
       var_list[1].setRange("signal", 103, 2286)
       ##plot the data
       roohist_pesudodata.plotOn(mjj_frame,
                                #ROOT.RooFit.Cut("Mj1>50 && Mj1<500 && Mj2>103 && Mj2<2286"),
                                #ROOT.RooFit.Range("signal")
                                )    
       ## plot the total pdf
       model_b.plotOn(mjj_frame, 
                      #ROOT.RooFit.ProjWData(RooArgSet(var_list[0], var_list[1]),roohist_pesudodata),  
                      #ROOT.RooFit.Range("signal"),
                      #RooFit.ProjectionRange("signal"),
                      RooFit.LineColor(ROOT.kBlack),
                      RooFit.LineWidth(1)
                      )
       ##plot the divid hist
       hpull = mjj_frame.pullHist() 
       mjj_Pull_Frame = var_list[i].frame()
       mjj_Pull_Frame.addPlotable(hpull,"P")

       ## plot the components pdf
       for sample in bkg_samples:
         model_b.plotOn(mjj_frame,ROOT.RooFit.Components(pdf_names[sample][category]["nominal"].GetName()),ROOT.RooFit.LineColor(samples_color[sample]),ROOT.RooFit.LineWidth(1))
       #model_b.plotOn(mjj_frame,
                      #ROOT.RooFit.ProjectionRange("signal"),
       #               ROOT.RooFit.Components(PDF_list[4].GetName()),ROOT.RooFit.LineColor(ROOT.kViolet),ROOT.RooFit.LineWidth(1))
       #roohist_signal.plotOn(mjj_frame,ROOT.RooFit.Cut("Mj1>50 && Mj1<500 && Mj2>103 && Mj2<2286"), ROOT.RooFit.DrawOption("P"), ROOT.RooFit.LineWidth(1),ROOT.RooFit.LineColor(13),ROOT.RooFit.MarkerColor(13), ROOT.RooFit.LineStyle(2))

       canvas.cd(1).SetPad(0.01, 0.23, 0.99, 0.98)
       gPad.SetLeftMargin(0.1) 
       gPad.SetLogy()

       # Get the maxmum of y axis
       max_bin   = hist1D_names["JetHT"][category]["nominal"]["mjj"].GetMaximumBin()
       max_value = hist1D_names["JetHT"][category]["nominal"]["mjj"].GetBinContent(max_bin)
       #print (max_value)

       # Draw the frame
       mjj_frame.SetTitle("")
       mjj_frame.GetXaxis().SetTitleOffset(0.9)
       mjj_frame.GetXaxis().SetTitleSize(0.050)
       mjj_frame.GetXaxis().SetLabelSize(0)
       #frame.GetYaxis().SetTitle("Events")
       mjj_frame.GetYaxis().SetLabelSize(0.050)
       mjj_frame.GetYaxis().SetTitleSize(0.050)
       mjj_frame.GetYaxis().SetTitleOffset(0.9)
       mjj_frame.SetMaximum(max_value*10)
       mjj_frame.SetMinimum(1)
       mjj_frame.Draw()
   
       canvas.cd(2).SetPad(0.01, 0.03, 0.99, 0.29)
       gPad.SetLeftMargin(0.1)
       gPad.SetBottomMargin(0.27)
   
       mjj_Pull_Frame.SetTitle("")
       mjj_Pull_Frame.GetXaxis().SetTitle("Mvv [GeV]")
       mjj_Pull_Frame.GetXaxis().SetTitleOffset(0.9)
       mjj_Pull_Frame.GetXaxis().SetTitleSize(0.15)
       mjj_Pull_Frame.GetXaxis().SetLabelSize(0.15)
       #frame.GetYaxis().SetTitle("Events")
       mjj_Pull_Frame.GetYaxis().SetTitle("Pull")
       mjj_Pull_Frame.GetYaxis().SetLabelSize(0.150)
       mjj_Pull_Frame.GetYaxis().SetTitleSize(0.150)
       mjj_Pull_Frame.GetYaxis().SetTitleOffset(0.9)
       mjj_Pull_Frame.Draw()

       line = ROOT.TLine(min_mjj, 0, max_mjj, 0)
       line.SetLineColor(2)
       line.Draw("same")
   
     canvas.cd(1)
     pad = ROOT.TPad("pad","pad",0.01,0.01,0.99,0.99)
     ROOT.gPad.RedrawAxis()
     channelText = ""
     channelTextFont   = 42
     channelTextSize   = 0.060
     cmsText     = "CMS"
     cmsTextFont   = 61  # default is helvetic-bold
     #writeExtraText = true
     extraText   = "Simulation"
     extraTextFont = 52  # default is helvetica-italics
     #text sizes and text offsets with respect to the top frame in unit of the top margin size
     lumiTextSize     = 0.7
     lumiTextOffset   = 0.2
     cmsTextSize      = 0.75
     cmsTextOffset    = 0.1  # only used in outOfFrame version
     relPosX    = 0.045
     relPosY    = 0.035
     relExtraDY = 1.2
     # ratio of "CMS" and extra text size
     extraOverCmsTextSize  = 0.65
     lumi_13TeV = "41.5 fb^{-1}"
     #lumiText += lumi_13TeV
     lumiText ="2017"+" (13 TeV)"
     t = pad.GetTopMargin()
     b = pad.GetBottomMargin()
     r = pad.GetRightMargin()
     l = pad.GetLeftMargin()
     latex = ROOT.TLatex()
     latex.SetNDC()
     latex.SetTextAngle(0)
     latex.SetTextColor(ROOT.kBlack)
     extraTextSize = extraOverCmsTextSize*cmsTextSize
     latex.SetTextFont(42)
     latex.SetTextAlign(31)
     latex.SetTextSize(lumiTextSize*t)
     latex.DrawLatex(1-r,0.92,lumiText)
     latex.SetTextFont(cmsTextFont)
     latex.SetTextAlign(11)
     latex.SetTextSize(cmsTextSize*t)
     latex.DrawLatex(l+0.01, 0.92,cmsText)
     #latex.DrawLatex(l+0.045, 1-t+lumiTextOffset*t-0.08,cmsText)
     latex.SetTextFont(extraTextFont)
     latex.SetTextSize(extraTextSize*t)  
     latex.DrawLatex(l+0.12, 0.92, extraText)  
  
     legend = ROOT.TLegend(0.65, 0.55, 0.85, 0.89)
     legend.SetTextSize(0.040)
     legend.SetTextFont(62)
     legend.SetFillColor(0)
     legend.SetBorderSize(0)
     legend.AddEntry(hist1D_names["JetHT"][category]["nominal"]["mjj"], "Pseudo data", "PE")
     for sample in bkg_samples:
       legend.AddEntry(hist1D_names[sample][category]["nominal"]["mjj"], sample, "l")
     legend.Draw()

     canvas.Draw()
     if i == 0:
       canvas.SaveAs("Post_fit_b_mj1.pdf")
     if i == 1:
       canvas.SaveAs("Post_fit_b_mj2.pdf")
     if i == 2:
       canvas.SaveAs("Post_fit_b_mjj.pdf")

     del pad


def plot_prefit(Cut,Variable,bkg_samples,samples_color,categories,varname_list,hist_names):
   for category in categories:
     #print (sample)
     print (category)
     old_canvas = ROOT.gROOT.FindObject("canvas")
     if old_canvas:
        old_canvas.Close()  
     canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
     hStack = ROOT.THStack()
     c1_1 = ROOT.TPad("c1_1", "newpad", 0.01,0.01,0.99,0.32)
     c1_1.Draw();
     c1_1.cd();
     c1_1.SetTopMargin(0.045);
     c1_1.SetBottomMargin(0.3);
     c1_1.SetRightMargin(0.035);
     c1_1.SetLeftMargin(0.11);
     
     gDATA = ROOT.TGraphAsymmErrors(hist_names["JetHT"][category]["nominal"])
     gDATA.SetLineWidth(2)
     gDATA.SetLineColor(1)
     gDATA.SetMarkerColor(1)
     gDATA.SetMarkerStyle(20)
     gDATA.SetMarkerSize(1.3)

     hist_names["JetHT"][category]["nominal"].SetLineWidth(2)
     hist_names["JetHT"][category]["nominal"].SetLineColor(1)
     hist_names["JetHT"][category]["nominal"].SetMarkerColor(1)
     hist_names["JetHT"][category]["nominal"].SetMarkerStyle(20)
     hist_names["JetHT"][category]["nominal"].SetMarkerSize(1.3)

     RATIO2 = hist_names["JetHT"][category]["nominal"].Clone("RATIO2") 
     for i, sample in enumerate(bkg_samples):
       if sample == "XToYYprime_MX3000":
         continue
       print (sample)
       hist_names[sample][category]["nominal"].SetFillColor(samples_color[sample])
       hist_names[sample][category]["nominal"].SetLineColor(samples_color[sample])
       if i == 0:
         hist_bkg = hist_names[sample][category]["nominal"].Clone("hist_total")
         hStack.Add(hist_names[sample][category]["nominal"])
       else:
         hist_bkg.Add(hist_names[sample][category]["nominal"])
         hStack.Add(hist_names[sample][category]["nominal"])
     
     n_bins = hist_names[sample][category]["nominal"].GetNbinsX()
     x = array('d', [])
     y = array('d', [])
     #y_up = array('d', [])
     #y_down = array('d', [])
     ex = array('d', [])
     ey_high = array('d', [])
     ey_low = array('d', [])

     for i in range(1, n_bins + 1):
       x.append(hist_bkg.GetBinCenter(i))
       ex.append(hist_bkg.GetBinWidth(i) / 2)
       RATIO2.SetBinContent(i,1)
       if (hist_bkg.GetBinContent(i)!=0):
         y.append(hist_names["JetHT"][category]["nominal"].GetBinContent(i)/hist_bkg.GetBinContent(i))    
         ey_low.append(math.sqrt(gDATA.GetErrorYlow(i)*gDATA.GetErrorYlow(i))/hist_bkg.GetBinContent(i))
         ey_high.append(math.sqrt(gDATA.GetErrorYhigh(i)*gDATA.GetErrorYhigh(i))/hist_bkg.GetBinContent(i))
         RATIO2.SetBinError(i,hist_bkg.GetBinError(i)/hist_bkg.GetBinContent(i))
       else:
         y.append(-1)
         ey_low.append(0)
         ey_high.append(0)

     min_x = hist_bkg.GetXaxis().GetXmin()
     max_x = hist_bkg.GetXaxis().GetXmax()
      
     RATIO = ROOT.TGraphAsymmErrors(n_bins, x, y, ex, ex, ey_low, ey_high)
     RATIO.Draw("AE0p")
     RATIO.SetMarkerColor(1)
     RATIO.SetMarkerStyle(21)
     RATIO.SetMarkerSize(1.0)
     RATIO.SetMaximum(1.3)
     RATIO.SetMinimum(0.7)
     RATIO.SetLineColor(1)
     RATIO.SetLineWidth(2)
     RATIO.SetTitle("")

     RATIO.GetYaxis().SetTitle("observed/expected")
     RATIO.GetYaxis().SetNdivisions(505)
     RATIO.GetYaxis().SetTitleSize(0.09)
     RATIO.GetYaxis().SetTitleOffset(0.5)
     RATIO.GetYaxis().SetLabelSize(0.13)

     RATIO.GetXaxis().SetTitle(Variable+" [GeV]")
     RATIO.GetXaxis().SetTitleSize(0.15)
     RATIO.GetXaxis().SetTitleOffset(0.9)
     RATIO.GetXaxis().SetLabelSize(0.13)

     RATIO.GetXaxis().SetRangeUser(min_x,max_x)
     RATIO2.SetFillStyle(3002)
     RATIO2.SetFillColor(12)
     RATIO2.SetLineColor(12)
     RATIO2.SetMarkerSize(0)
     RATIO2.Draw("E2same")
     RATIO.Draw("E0psame")

     line = ROOT.TLine(min_x, 1, max_x, 1)
     line.SetLineColor(1)
     line.Draw()

     canvas.cd()
     c1_2 = ROOT.TPad("c1_2", "newpad", 0.01,0.32,0.99,0.99)
     c1_2.Draw()
     c1_2.cd()
     c1_2.SetTopMargin(0.08)
     c1_2.SetBottomMargin(0.02)
     c1_2.SetRightMargin(0.035)
     c1_2.SetLeftMargin(0.11)
     c1_2.SetLogy()
     
     hStack.Draw("histo")
     hStack.SetMinimum(10)
     hStack.GetYaxis().SetTitleSize(0.070)
     hStack.GetXaxis().SetTitleSize(0.070)
     hStack.GetYaxis().SetLabelSize(0.070)
     hStack.GetXaxis().SetLabelSize(0.0)
     hStack.SetTitle("")
     hStack.GetYaxis().SetTitle("Events")
     hStack.GetXaxis().SetTitle("")
     hStack.GetYaxis().SetTitleOffset(0.80)
     hStack.GetXaxis().SetTitleOffset(0.85)
     
     hist_names["JetHT"][category]["nominal"].Draw("E same");
     
     hist_bkg.SetFillStyle(3005)  
     hist_bkg.SetFillColor(12)   
     hist_bkg.SetLineColor(12)    
     hist_bkg.Draw("E2same")

     max_bin = hist_bkg.GetMaximumBin()
     max_value = hist_bkg.GetBinContent(max_bin)

     hStack.SetMaximum(max_value*10)

     legend = ROOT.TLegend(0.70, 0.55, 0.89, 0.89)
     legend.SetTextSize(0.050)
     legend.SetTextFont(62)
     legend.SetFillColor(0)
     for i, sample in enumerate(bkg_samples[::-1]):
       legend.AddEntry(hist_names[sample][category]["nominal"],sample, "F")
     legend.SetBorderSize(0)
     legend.Draw()

     canvas.cd()
     pad = ROOT.TPad("pad","pad",0.01,0.01,0.99,0.99)
     ROOT.gPad.RedrawAxis()
     channelText = ""
     channelTextFont   = 42
     channelTextSize   = 0.060
     cmsText     = "CMS"
     cmsTextFont   = 61  # default is helvetic-bold
     #writeExtraText = true
     extraText   = "Simulation"
     extraTextFont = 52  # default is helvetica-italics
     #text sizes and text offsets with respect to the top frame in unit of the top margin size
     lumiTextSize     = 0.5
     lumiTextOffset   = 0.2
     cmsTextSize      = 0.55
     cmsTextOffset    = 0.1  # only used in outOfFrame version
     relPosX    = 0.045
     relPosY    = 0.035
     relExtraDY = 1.2
     # ratio of "CMS" and extra text size
     extraOverCmsTextSize  = 0.65
     lumi_13TeV = "41.5 fb^{-1}"
     #lumiText += lumi_13TeV
     lumiText ="2017"+" (13 TeV)"
     t = pad.GetTopMargin()
     b = pad.GetBottomMargin()
     r = pad.GetRightMargin()
     l = pad.GetLeftMargin()
     latex = ROOT.TLatex()
     latex.SetNDC()
     latex.SetTextAngle(0)
     latex.SetTextColor(ROOT.kBlack)
     extraTextSize = extraOverCmsTextSize*cmsTextSize
     latex.SetTextFont(42)
     latex.SetTextAlign(31)
     latex.SetTextSize(lumiTextSize*t)
     latex.DrawLatex(1-r,0.95,lumiText)
     latex.SetTextFont(cmsTextFont)
     latex.SetTextAlign(11)
     latex.SetTextSize(cmsTextSize*t)
     latex.DrawLatex(l+0.01, 0.95,cmsText)
     #latex.DrawLatex(l+0.045, 1-t+lumiTextOffset*t-0.08,cmsText)
     latex.SetTextFont(extraTextFont)
     latex.SetTextSize(extraTextSize*t)  
     latex.DrawLatex(l+0.12, 0.95, extraText)  

     canvas.Draw()
     canvas.SaveAs("plots_prefit/"+Cut+"_"+Variable+"_"+category+"_prefit.pdf")
   ROOT.gDirectory.GetList().Remove(canvas)             


def Add_generator_shower_sys(categories,hist_covert3Dto1D,hist_3D,hist_1D):
   ##When you plot the sys, you should use the full 3D and 1D bins, but when do the combine just a part of bins, the 1D bins should be just thousands.
   for category in categories:  
     n_bins = hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["nominal"].GetNbinsX()
     ###shower model sys
     hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["showerUp"]   = hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8"+"_"+category+"_showerUp")
     hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["showerDown"] = hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8"+"_"+category+"_showerDown")
     ###ME model sys
     hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEUp"]   = hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8"+"_"+category+"_MEUp")
     hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEDown"] = hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8"+"_"+category+"_MEDown")
     ###ME and shower model sys
     hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEshowerUp"]   = hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8"+"_"+category+"_MEshowerUp")
     hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEshowerDown"] = hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8"+"_"+category+"_MEshowerDown")

     for i_bin in range(1, n_bins+1):
       content_shower   = hist_covert3Dto1D["QCD_madgraph_herwig7"][category]["nominal"].GetBinContent(i_bin) 
       content_ME       = hist_covert3Dto1D["QCD_pythia8_Pt"][category]["nominal"].GetBinContent(i_bin) 
       content_MEshower = hist_covert3Dto1D["QCD_herwig7_Pt"][category]["nominal"].GetBinContent(i_bin) 

       # Make the Up sys, down sys is same
       content_showerUp   = content_shower
       content_MEUp       = content_ME
       content_MEshowerUp = content_MEshower

       hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["showerUp"].SetBinContent(i_bin,content_showerUp) 
       hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEUp"].SetBinContent(i_bin,content_MEUp) 
       hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["MEshowerUp"].SetBinContent(i_bin,content_MEshowerUp) 
       #print(f"sys shower up is {content_showerUp}")
       #print(f"sys ME up is {content_MEUp}")
       #print(f"sys MEshower up is {content_MEshowerUp}")
       #print(" ")
  
     ### convert the 1D to 3D
     ###shower model sys
     hist_3D["QCD_madgraph_pythia8"][category]["showerUp"]   = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8_3D"+"_"+category+"_showerUp")
     hist_3D["QCD_madgraph_pythia8"][category]["showerDown"] = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8_3D"+"_"+category+"_showerDown")
     ###ME model sys
     hist_3D["QCD_madgraph_pythia8"][category]["MEUp"]   = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8_3D"+"_"+category+"_MEUp")
     hist_3D["QCD_madgraph_pythia8"][category]["MEDown"] = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8_3D"+"_"+category+"_MEDown")
     ###ME and shower model sys
     hist_3D["QCD_madgraph_pythia8"][category]["MEshowerUp"]   = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8_3D"+"_"+category+"_MEshowerUp")
     hist_3D["QCD_madgraph_pythia8"][category]["MEshowerDown"] = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].Clone("QCD_madgraph_pythia8_3D"+"_"+category+"_MEshowerDown")
    

     projX = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].ProjectionX("projX")
     projY = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].ProjectionY("projY") 
     projZ = hist_3D["QCD_madgraph_pythia8"][category]["nominal"].ProjectionZ("projZ")
     nBinsX = projX.GetNbinsX()
     nBinsY = projY.GetNbinsX()  
     nBinsZ = projZ.GetNbinsX()

     #print (f"hist 3D number of x,y,z bins is {nBinsX},{nBinsY},{nBinsZ}")

     # 遍历所有可能的(i,j,k)组合
     for i in range(1, nBinsX + 1):
         for j in range(1, nBinsY + 1):
             for k in range(1, nBinsZ + 1):
                 # 计算对应的一维索引（与原始转换公式相同）
                 index_1d = k + nBinsZ * ((j-1) + nBinsY * (i-1))
                 showerUp_content   = hist_3D["QCD_madgraph_herwig7"][category]["nominal"].GetBinContent(i,j,k)
                 MEUp_content       = hist_3D["QCD_pythia8_Pt"][category]["nominal"].GetBinContent(i,j,k)
                 MEshowerUp_content = hist_3D["QCD_herwig7_Pt"][category]["nominal"].GetBinContent(i,j,k)
                 hist_3D["QCD_madgraph_pythia8"][category]["showerUp"].SetBinContent(i, j, k, showerUp_content)
                 hist_3D["QCD_madgraph_pythia8"][category]["MEUp"].SetBinContent(i, j, k, MEUp_content)
                 hist_3D["QCD_madgraph_pythia8"][category]["MEshowerUp"].SetBinContent(i, j, k, MEshowerUp_content)
     

     ###shower model sys
     hist_1D["QCD_madgraph_pythia8"][category]["showerUp"]["fatjet"] = hist_3D["QCD_madgraph_pythia8"][category]["showerUp"].ProjectionX("hist_QCD_madgraph_pythia8"+"_"+category+"_showerUp_fatjet")
     hist_1D["QCD_madgraph_pythia8"][category]["showerUp"]["2jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["showerUp"].ProjectionY("hist_QCD_madgraph_pythia8"+"_"+category+"_showerUp_2jets" )
     hist_1D["QCD_madgraph_pythia8"][category]["showerUp"]["3jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["showerUp"].ProjectionZ("hist_QCD_madgraph_pythia8"+"_"+category+"_showerUp_3jets" )
                                                                                                            
     hist_1D["QCD_madgraph_pythia8"][category]["showerDown"]["fatjet"] = hist_3D["QCD_madgraph_pythia8"][category]["showerDown"].ProjectionX("hist_QCD_madgraph_pythia8"+"_"+category+"_showerDown_fatjet")
     hist_1D["QCD_madgraph_pythia8"][category]["showerDown"]["2jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["showerDown"].ProjectionY("hist_QCD_madgraph_pythia8"+"_"+category+"_showerDown_2jets" )
     hist_1D["QCD_madgraph_pythia8"][category]["showerDown"]["3jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["showerDown"].ProjectionZ("hist_QCD_madgraph_pythia8"+"_"+category+"_showerDown_3jets" )
     ###ME model sys
     hist_1D["QCD_madgraph_pythia8"][category]["MEUp"]["fatjet"] = hist_3D["QCD_madgraph_pythia8"][category]["MEUp"].ProjectionX("hist_QCD_madgraph_pythia8"+"_"+category+"_MEUp_fatjet")
     hist_1D["QCD_madgraph_pythia8"][category]["MEUp"]["2jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["MEUp"].ProjectionY("hist_QCD_madgraph_pythia8"+"_"+category+"_MEUp_2jets" )
     hist_1D["QCD_madgraph_pythia8"][category]["MEUp"]["3jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["MEUp"].ProjectionZ("hist_QCD_madgraph_pythia8"+"_"+category+"_MEUp_3jets" )

     hist_1D["QCD_madgraph_pythia8"][category]["MEDown"]["fatjet"] = hist_3D["QCD_madgraph_pythia8"][category]["MEDown"].ProjectionX("hist_QCD_madgraph_pythia8"+"_"+category+"_MEDown_fatjet")
     hist_1D["QCD_madgraph_pythia8"][category]["MEDown"]["2jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["MEDown"].ProjectionY("hist_QCD_madgraph_pythia8"+"_"+category+"_MEDown_2jets" )
     hist_1D["QCD_madgraph_pythia8"][category]["MEDown"]["3jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["MEDown"].ProjectionZ("hist_QCD_madgraph_pythia8"+"_"+category+"_MEDown_3jets" )

     ###ME and shower model sys
     hist_1D["QCD_madgraph_pythia8"][category]["MEshowerUp"]["fatjet"] = hist_3D["QCD_madgraph_pythia8"][category]["MEshowerUp"].ProjectionX("hist_QCD_madgraph_pythia8"+"_"+category+"_MEshowerUp_fatjet")
     hist_1D["QCD_madgraph_pythia8"][category]["MEshowerUp"]["2jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["MEshowerUp"].ProjectionY("hist_QCD_madgraph_pythia8"+"_"+category+"_MEshowerUp_2jets" )
     hist_1D["QCD_madgraph_pythia8"][category]["MEshowerUp"]["3jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["MEshowerUp"].ProjectionZ("hist_QCD_madgraph_pythia8"+"_"+category+"_MEshowerUp_3jets" )

     hist_1D["QCD_madgraph_pythia8"][category]["MEshowerDown"]["fatjet"] = hist_3D["QCD_madgraph_pythia8"][category]["MEshowerDown"].ProjectionX("hist_QCD_madgraph_pythia8"+"_"+category+"_MEshowerDown_fatjet")
     hist_1D["QCD_madgraph_pythia8"][category]["MEshowerDown"]["2jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["MEshowerDown"].ProjectionY("hist_QCD_madgraph_pythia8"+"_"+category+"_MEshowerDown_2jets" )
     hist_1D["QCD_madgraph_pythia8"][category]["MEshowerDown"]["3jets" ] = hist_3D["QCD_madgraph_pythia8"][category]["MEshowerDown"].ProjectionZ("hist_QCD_madgraph_pythia8"+"_"+category+"_MEshowerDown_3jets" )

     #print(f"sys madgraph pythia8 is ", hist_covert3Dto1D["QCD_madgraph_pythia8"][category]["nominal"].Integral())
     #print(f"sys madgraph herwig7 is ", hist_covert3Dto1D["QCD_madgraph_herwig7"][category]["nominal"].Integral())
     #print(f"sys herwig7 is ", hist_covert3Dto1D["QCD_herwig7_Pt"][category]["nominal"].Integral())
     #print(f"sys pythia8 is ", hist_covert3Dto1D["QCD_pythia8_Pt"][category]["nominal"].Integral())
     #print(" ")
     

def functional_sys_weight(dimension,sys_type,mass):
   weight = 1
   if dimension == "fatjet":
     if sys_type == "up":
       weight = 1+(mass-200)*0.2/500
     if sys_type == "down":
       weight = 1-(mass-200)*0.2/500
     if sys_type == "nominal":
       weight = 1
   if dimension == "2jets":
     if sys_type == "up":
       weight = 1+(mass-1000)*0.2/3300
     if sys_type == "down":
       weight = 1-(mass-1000)*0.2/3300
     if sys_type == "nominal":
       weight = 1
   if dimension == "3jets":
     if sys_type == "up":
       weight = 1+(mass-2050)*0.4/5200
     if sys_type == "down":
       weight = 1-(mass-2050)*0.4/5200
     if sys_type == "nominal":
       weight = 1
   return weight

def functional_inv_sys_weight(dimension,sys_type,mass):
   weight = 1
   if dimension == "fatjet":
     if sys_type == "up":
       weight = 1.2-(30/mass)
     if sys_type == "down":
       weight = 0.8+(30/mass)
     if sys_type == "nominal":
       weight = 1
   if dimension == "2jets":
     if sys_type == "up":
       weight = 1.2-(150*0.2/mass)
     if sys_type == "down":
       weight = 0.8+(150*0.2/mass)
     if sys_type == "nominal":
       weight = 1
   if dimension == "3jets":
     if sys_type == "up":
       weight = 1.2-(405/mass)
     if sys_type == "down":
       weight = 0.8+(405/mass)
     if sys_type == "nominal":
       weight = 1
   return weight

def Add_functional_sys(sample,category,hist3D_names,hist1D_names,varname_list):
  h3 = hist3D_names[sample][category]["nominal"] 
  for vn in varname_list:
    for dn in ["Up","Down","invUp","invDown"]:
      hist3D_names[sample][category]["m"+vn+dn]   = hist3D_names[sample][category]["nominal"].Clone(sample+"_3D_"+category+"_m"+vn+dn)
  
  projX = h3.ProjectionX(sample+"_functionalsys_projX")
  projY = h3.ProjectionY(sample+"_functionalsys_projY") 
  projZ = h3.ProjectionZ(sample+"_functionalsys_projZ")

  nBinsX = projX.GetNbinsX()
  nBinsY = projY.GetNbinsX()  
  nBinsZ = projZ.GetNbinsX()

  # 遍历所有可能的(i,j,k)组合,fill in the three dimension
  for i in range(1, nBinsX + 1):
      for j in range(1, nBinsY + 1):
          for k in range(1, nBinsZ + 1):
              masses= [ projX.GetXaxis().GetBinUpEdge(i), projY.GetXaxis().GetBinUpEdge(j), projZ.GetXaxis().GetBinUpEdge(k)]
              nominal_value = hist3D_names[sample][category]["nominal"].GetBinContent(i,j,k)
              for vn, mass in zip(varname_list,masses):
                hist3D_names[sample][category]["m"+vn+"Up"].SetBinContent(i,j,k,nominal_value*functional_sys_weight(vn,"up",mass))
                hist3D_names[sample][category]["m"+vn+"Down"].SetBinContent(i,j,k,nominal_value*functional_sys_weight(vn,"down",mass))
                hist3D_names[sample][category]["m"+vn+"invUp"].SetBinContent(i,j,k,nominal_value*functional_inv_sys_weight(vn,"up",mass))
                hist3D_names[sample][category]["m"+vn+"invDown"].SetBinContent(i,j,k,nominal_value*functional_inv_sys_weight(vn,"down",mass))
   

  ### functional sys up and down
  for vn in varname_list:
   for dn in ["Up","Down","invUp","invDown"]:
     hist1D_names[sample][category]["m"+vn+dn]["fatjet"] = hist3D_names[sample][category]["m"+vn+dn].ProjectionX("hist_"+sample+"_"+category+"_m"+vn+dn+"_mfatjet")
     hist1D_names[sample][category]["m"+vn+dn]["2jets"] = hist3D_names[sample][category]["m"+vn+dn].ProjectionY("hist_"+sample+"_"+category+"_m"+vn+dn+"_m2jets")
     hist1D_names[sample][category]["m"+vn+dn]["3jets"] = hist3D_names[sample][category]["m"+vn+dn].ProjectionZ("hist_"+sample+"_"+category+"_m"+vn+dn+"_m3jets")

def plot_QCD_diff_generator_shower(categories,varname_list,hist1D_names):
   ##When you plot the sys, you should use the full 3D and 1D bins, but when do the combine just a part of bins, the 1D bins should be just thousands.
   for category in categories:
     for var in varname_list:
       old_canvas = ROOT.gROOT.FindObject("canvas")
       if old_canvas:
         old_canvas.Close()  
       canvas = ROOT.TCanvas("canvas", "Stacked PDF Projection", 800, 600)
       #hist1D_names["QCD_madgraph_herwig7"][category]["nominal"][var].Scale(1500)

       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].SetLineColor(ROOT.kBlack)
       hist1D_names["QCD_herwig7_Pt"][category]["nominal"][var].SetLineColor(ROOT.kRed)
       hist1D_names["QCD_pythia8_Pt"][category]["nominal"][var].SetLineColor(ROOT.kBlue)
       hist1D_names["QCD_madgraph_herwig7"][category]["nominal"][var].SetLineColor(ROOT.kGreen)
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].SetFillColor(0)
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].SetTitle("")
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].GetYaxis().SetTitleSize(0.045);
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].GetXaxis().SetTitleSize(0.045);
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].GetYaxis().SetLabelSize(0.045);
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].GetXaxis().SetLabelSize(0.045);
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].GetYaxis().SetTitleOffset(1.2);
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].GetYaxis().SetTitle("Number of Events");
       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].GetXaxis().SetTitle(var+" mass [GeV]");

       hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var].Draw("HIST")
       hist1D_names["QCD_madgraph_herwig7"][category]["nominal"][var].Draw("HIST same")
       hist1D_names["QCD_herwig7_Pt"][category]["nominal"][var].Draw("HIST same")
       hist1D_names["QCD_pythia8_Pt"][category]["nominal"][var].Draw("HIST same")

       legend = ROOT.TLegend(0.50, 0.75, 0.88, 0.89)
       legend.SetTextSize(0.040)
       legend.SetTextFont(62)
       legend.SetFillColor(0)
       legend.AddEntry(hist1D_names["QCD_madgraph_pythia8"][category]["nominal"][var],"QCD madgraph pythia8", "l")
       legend.AddEntry(hist1D_names["QCD_madgraph_herwig7"][category]["nominal"][var],"QCD madgraph herwig7", "l")
       legend.AddEntry(hist1D_names["QCD_herwig7_Pt"][category]["nominal"][var],"QCD herwig7", "l")
       legend.AddEntry(hist1D_names["QCD_pythia8_Pt"][category]["nominal"][var],"QCD pythia8", "l")
       legend.SetBorderSize(0)
       legend.Draw()

       canvas.Draw()
       canvas.SaveAs("plots_sys/QCD_madgraph_pythia8_"+category+"_ME_and_shower_"+var+"_sys.pdf")
   ROOT.gDirectory.GetList().Remove(canvas)             
    


def plot_sys(store_path,samples,categories,systematics_names,varname_list,hist1D_names):
   CreatDirectory(store_path+"/plots_sys")
   for sample in samples:
     for category in categories:
       for systematic in systematics_names:
         for var in varname_list:
           if hist1D_names[sample][category][systematic+"Down"][var] == None:
              continue
           print (sample)
           print (category)
           old_canvas = ROOT.gROOT.FindObject("canvas_sys")
           if old_canvas:
              old_canvas.Close()  
           canvas_sys = ROOT.TCanvas("canvas_sys", "Stacked PDF Projection", 800, 600)
           c1_1 = ROOT.TPad("c1_1", "newpad", 0.01,0.01,0.99,0.32)
           c1_1.Draw();
           c1_1.cd();
           c1_1.SetTopMargin(0.045);
           c1_1.SetBottomMargin(0.3);
           c1_1.SetRightMargin(0.035);
           c1_1.SetLeftMargin(0.11);

           n_bins = hist1D_names[sample][category]["nominal"][var].GetNbinsX()
           x = array('d', [])
           y_up = array('d', [])
           y_down = array('d', [])
           ex = array('d', [])
           ey_high = array('d', [])
           ey_low = array('d', [])

           min_x = hist1D_names[sample][category]["nominal"][var].GetXaxis().GetXmin()
           max_x = hist1D_names[sample][category]["nominal"][var].GetXaxis().GetXmax()

           print (min_x, max_x)
           for i in range(1, n_bins + 1):
             nom = hist1D_names[sample][category]["nominal"][var].GetBinContent(i)
             up = hist1D_names[sample][category][systematic+"Up"][var].GetBinContent(i)
             down = hist1D_names[sample][category][systematic+"Down"][var].GetBinContent(i)
             x.append(hist1D_names[sample][category]["nominal"][var].GetBinCenter(i))
             ex.append(hist1D_names[sample][category]["nominal"][var].GetBinWidth(i) / 2)
 
             if nom != 0:
                 y_up.append(up/nom)  
                 y_down.append(down/nom)  
                 ey_high.append(0.0)
                 ey_low.append(0.0)
             else:
                 y_up.append(0.0)
                 y_down.append(0.0)
                 ey_high.append(0.0)
                 ey_low.append(0.0)

           ratio_up = ROOT.TGraphAsymmErrors(n_bins, x, y_up, ex, ex, ey_low, ey_high)
           ratio_up.Draw("AE0p")
           ratio_up.SetMarkerColor(1)
           ratio_up.SetMarkerStyle(21)
           ratio_up.SetMarkerSize(0.0)
           ratio_up.SetMaximum(1.2)
           ratio_up.SetMinimum(0.8)
           ratio_up.SetLineColor(2)
           ratio_up.SetLineWidth(2)
           ratio_up.GetXaxis().SetTitleOffset(0.9)
           ratio_up.GetYaxis().SetTitleOffset(0.5)
           ratio_up.SetTitle("")
           ratio_up.GetYaxis().SetTitle("ratio")
           ratio_up.GetXaxis().SetTitle(var+" mass [GeV]")
           ratio_up.GetXaxis().SetLabelSize(0.13)
           ratio_up.GetYaxis().SetLabelSize(0.13)
           ratio_up.GetXaxis().SetTitleSize(0.15)
           ratio_up.GetYaxis().SetTitleSize(0.09)
           ratio_up.GetYaxis().SetNdivisions(505)
           ratio_up.GetXaxis().SetRangeUser(min_x,max_x)

           ratio_down = ROOT.TGraphAsymmErrors(n_bins, x, y_down, ex, ex, ey_low, ey_high)
           ratio_down.SetMarkerColor(1)
           ratio_down.SetMarkerStyle(21)
           ratio_down.SetMarkerSize(0.0)
           #ratio_down.SetMaximum(1.2)
           #ratio_down.SetMinimum(0.8)
           ratio_down.SetLineColor(4)
           ratio_down.SetLineWidth(2)
           ratio_down.Draw("E0psame")
           ratio_up.Draw("E0psame")
           line = ROOT.TLine(hist1D_names[sample][category]["nominal"][var].GetXaxis().GetXmin(), 1, hist1D_names[sample][category]["nominal"][var].GetXaxis().GetXmax(), 1)
           line.SetLineColor(1)
           line.Draw()

           canvas_sys.cd()
           c1_2 = ROOT.TPad("c1_2", "newpad", 0.01,0.32,0.99,0.99)
           c1_2.Draw()
           c1_2.cd()
           c1_2.SetTopMargin(0.08);
           c1_2.SetBottomMargin(0.02);
           c1_2.SetRightMargin(0.035);
           c1_2.SetLeftMargin(0.11);
           hist1D_names[sample][category]["nominal"][var].SetLineColor(ROOT.kBlack)
           hist1D_names[sample][category]["nominal"][var].SetFillColor(0)
           hist1D_names[sample][category]["nominal"][var].SetTitle("")
           hist1D_names[sample][category]["nominal"][var].GetYaxis().SetTitleSize(0.045);
           hist1D_names[sample][category]["nominal"][var].GetXaxis().SetTitleSize(0.045);
           hist1D_names[sample][category]["nominal"][var].GetYaxis().SetLabelSize(0.045);
           hist1D_names[sample][category]["nominal"][var].GetXaxis().SetLabelSize(0.0);
           hist1D_names[sample][category]["nominal"][var].GetYaxis().SetTitleOffset(1.2);
           hist1D_names[sample][category]["nominal"][var].GetYaxis().SetTitle("Number of Events");

           hist1D_names[sample][category][systematic+"Up"][var].SetLineColor(ROOT.kRed)
           hist1D_names[sample][category][systematic+"Up"][var].SetFillColor(0)
           hist1D_names[sample][category][systematic+"Down"][var].SetLineColor(ROOT.kBlue)
           hist1D_names[sample][category][systematic+"Down"][var].SetFillColor(0)

           max_val = max(hist1D_names[sample][category][systematic+"Up"][var].GetMaximum(), hist1D_names[sample][category][systematic+"Down"][var].GetMaximum())
           hist1D_names[sample][category]["nominal"][var].SetMaximum(max_val*1.1)
           #max_bin = hist_bkg.GetMaximumBin()
           #max_value = hist_bkg.GetBinContent(max_bin)

           hist1D_names[sample][category]["nominal"][var].Draw("HIST")
           hist1D_names[sample][category][systematic+"Up"][var].Draw("HIST same")
           hist1D_names[sample][category][systematic+"Down"][var].Draw("HIST same")

           legend = ROOT.TLegend(0.50, 0.75, 0.88, 0.89)
           legend.SetTextSize(0.050)
           legend.SetTextFont(62)
           legend.SetFillColor(0)
           legend.AddEntry(hist1D_names[sample][category]["nominal"][var],sample+" nominal", "l")
           legend.AddEntry(hist1D_names[sample][category][systematic+"Up"][var],systematic+" up", "l")
           legend.AddEntry(hist1D_names[sample][category][systematic+"Down"][var],systematic+" down", "l")
           legend.SetBorderSize(0)
           legend.Draw()

           canvas_sys.Draw()
           canvas_sys.SaveAs(store_path+"/plots_sys/"+sample+"_"+category+"_"+systematic+"_"+var+"_mass.pdf")
   ROOT.gDirectory.GetList().Remove(canvas_sys)             

def MakePesudoData(signal_hist_list,bkg_hist_list):
   for i, signal_hist in enumerate(signal_hist_list):
     if i==0:
       hist_pesudo_data = signal_hist.Clone("hist_pesudo_data")
     else:
       hist_pesudo_data.Add(signal_hist)
   for j, bkg_hist in enumerate(bkg_hist_list):
     hist_pesudo_data.Add(bkg_hist) 

   return hist_pesudo_data

def MakePesudoData_bkgonly(bkg_hist_list):
   for j, bkg_hist in enumerate(bkg_hist_list):
     if j==0:
       hist_pesudo_data = bkg_hist.Clone("hist_pesudo_data")
     else:
       hist_pesudo_data.Add(bkg_hist)

   return hist_pesudo_data

def SumOtherBkg(otherbkg_hist_list):
   for i, bkg_hist in enumerate(otherbkg_hist_list):
     if i==0:
       hist_otherbkg = bkg_hist.Clone("hist_otherbkg")
     else:
       hist_otherbkg.Add(bkg_hist) 

   return hist_otherbkg

def MakeWorkspace_1D(category,samples,systematics,Roodatahist1D_names):
   w = ROOT.RooWorkspace("w")
   dataHist=Roodatahist1D_names["JetHT"][category]["nominal"]["mjj"] 
   print (category+ "data yield is ",  Roodatahist1D_names["JetHT"][category]["nominal"]["3jets"].sumEntries()) 

   getattr(w,'import')(dataHist)#,ROOT.RooFit.RenameVariable("data_obs","data_obs")
   for i, sample in enumerate(samples):  
     if sample in ["JetHT"]:
       continue
     for systematic in systematics:
       print (sample)
       if Roodatahist1D_names[sample][category][systematic]["mjj"] != None:
         getattr(w,'import')(Roodatahist1D_names[sample][category][systematic]["mjj"])#,ROOT.RooFit.RenameVariable(sample+"_"+category,sample+"_"+category)
   w.writeToFile("datacardInput_"+category+".root")
   
   print ("datacard " + category + " is done !")
   del w

def MakeWorkspace(category,signal_sample,bkg_samples,systematics,Roodatahist_names):
   w = ROOT.RooWorkspace("w")
   dataHist=Roodatahist_names["JetHT"][category]["nominal"] 

   getattr(w,'import')(dataHist,ROOT.RooFit.RenameVariable("data_obs","data_obs"))
   for i, bkg_sample in enumerate(bkg_samples):  
     if bkg_sample in ["JetHT"]:
       continue
     for systematic in systematics:
       print (bkg_sample)
       if Roodatahist_names[bkg_sample][category][systematic] != None:
         getattr(w,'import')(Roodatahist_names[bkg_sample][category][systematic])#,ROOT.RooFit.RenameVariable(sample+"_"+category,sample+"_"+category)
   for systematic in systematics:
     if Roodatahist_names[signal_sample][category][systematic] != None:
       getattr(w,'import')(Roodatahist_names[signal_sample][category][systematic])#,ROOT.RooFit.RenameVariable(sample+"_"+category,sample+"_"+category)
   w.writeToFile("datacardInput_"+category+".root")
   
   print ("datacard " + category + " is done !")
   del w

def MakeWorkspaceSimple(category,signal_sample,bkg_samples,systematics,datahist_names):
   w = ROOT.TFile.Open("datacardInputSimple_"+category+".root","RECREATE")
   dataHist=datahist_names["JetHT"][category]["nominal"] 

   w.WriteObject(dataHist,"data_obs")
   for i, bkg_sample in enumerate(bkg_samples):  
     if bkg_sample in ["JetHT"]:
       continue
     for systematic in systematics:
       print (bkg_sample)
       if datahist_names[bkg_sample][category][systematic] != None:
         w.WriteObject(datahist_names[bkg_sample][category][systematic],datahist_names[bkg_sample][category][systematic].GetName())
   for systematic in systematics:
     if datahist_names[signal_sample][category][systematic] != None:
       w.WriteObject(datahist_names[signal_sample][category][systematic],datahist_names[signal_sample][category][systematic].GetName())
   w.Close()
   
   print ("datacard " + category + " is done !")


def WriteDatacard(category,signal_sample,hist_names):
 mass_point = re.search(r'MX\d+', signal_sample).group()
 with open("datacard_"+category+"_"+mass_point+".txt", "w") as datacard:
  datacard.write("imax 1 number of bins\n") 
  datacard.write("jmax 6 number of backgrounds\n") 
  datacard.write("kmax * number of nuisance parameters\n") 
  datacard.write("------------------------------------\n") 
  #datacard.write("shapes * * datacardInput_%s.root w:$PROCESS w:$PROCESS_$SYSTEMATIC\n" % (category)) 
  datacard.write("shapes * * datacardInputSimple_%s.root $PROCESS $PROCESS_$SYSTEMATIC\n" % (category)) 
  datacard.write("------------------------------------\n") 
  datacard.write("bin bin1\n") 
  datacard.write("observation -1\n") 
  datacard.write("------------------------------------\n") 
  datacard.write("bin bin1 bin1 bin1 bin1 bin1 bin1 bin1\n") 
  datacard.write("process %s %s %s %s %s %s %s\n" % (hist_names[signal_sample][category]["nominal"].GetName(),hist_names["QCD_madgraph_pythia8"][category]["nominal"].GetName(),hist_names["TT"][category]["nominal"].GetName(),hist_names["ZJets"][category]["nominal"].GetName(),hist_names["WJets"][category]["nominal"].GetName(),hist_names["VV"][category]["nominal"].GetName(),hist_names["ST"][category]["nominal"].GetName()) )
  datacard.write("process 0 1 2 3 4 5 6\n")
  datacard.write("rate %f %f %f %f %f %f %f\n" % (hist_names[signal_sample][category]["nominal"].Integral(),hist_names["QCD_madgraph_pythia8"][category]["nominal"].Integral(),hist_names["TT"][category]["nominal"].Integral(),hist_names["ZJets"][category]["nominal"].Integral(),hist_names["WJets"][category]["nominal"].Integral(),hist_names["VV"][category]["nominal"].Integral(),hist_names["ST"][category]["nominal"].Integral())) 
  datacard.write("------------------------------------\n") 
  datacard.write("lumi lnN 1.025 - - - - - -\n") 
  datacard.write("QCD_norm lnN - 1.3 - - - - -\n") 
  datacard.write("ZJetsren_norm lnN - - - 1.01 - - -\n") 
  datacard.write("WJetsren_norm lnN - - - - 1.01 - -\n") 

  #datacard.write("isr_norm lnN - 1.01 1.05 1.05 1.05 1.05 1.01\n") 
  #datacard.write("fsr_norm lnN - 1.10 1.10 1.10 1.10 1.10 1.10\n") 

  datacard.write("QCD_madgraph_pythia8MEren  shape - 1 - - - - -\n") 
  datacard.write("TTMEren   shape - - 1 - - - -\n") 
  datacard.write("STMEren   shape - - - - - - 1\n") 

  datacard.write("QCD_madgraph_pythia8MEfac  shape - 1 - - - - -\n") 
  datacard.write("TTMEfac   shape - - 1 - - - -\n") 
  datacard.write("ZJetsMEfac shape - - - 1 - - -\n") 
  datacard.write("WJetsMEfac shape - - - - 1 - -\n") 
  datacard.write("STMEfac   shape - - - - - - 1\n") 

  datacard.write("shower  shape - 1 - - - - -\n") 
  datacard.write("ME  shape - 1 - - - - -\n") 
  datacard.write("MEshower  shape - 1 - - - - -\n") 

  datacard.write("mfatjet  shape - 1 - - - - -\n") 
  datacard.write("mfatjetinv  shape - 1 - - - - -\n") 
  datacard.write("m2jets  shape - 1 - - - - -\n") 
  datacard.write("m2jetsinv  shape - 1 - - - - -\n") 
  datacard.write("m3jets  shape - 1 - - - - -\n") 
  datacard.write("m3jetsinv  shape - 1 - - - - -\n") 

  datacard.write(signal_sample+"PSisr   shape 1 - - - - - -\n") 
  datacard.write("QCD_madgraph_pythia8PSisr   shape - 1 - - - - -\n") 
  datacard.write("TTPSisr   shape - - 1 - - - -\n") 
  datacard.write("ZJetsPSisr shape - - - 1 - - -\n") 
  datacard.write("WJetsPSisr shape - - - - 1 - -\n") 
  datacard.write("STPSisr   shape - - - - - - 1\n") 

  datacard.write("PSfsr   shape 1 1 1 1 1 1 1\n") 
  datacard.write("JuncTotal  shape 1 1 1 1 1 1 1\n") 

  #datacard.write("bin1 autoMCStats 0\n") 
