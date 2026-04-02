import ROOT
import os
import re
from array import array
import ctypes
from parameter import mj1_bins_3jets,mj2_bins_3jets,mjj_bins_3jets,mj1_bins_2fatjets,mj2_bins_2fatjets,mjj_bins_2fatjets
from parameter import signal_samples 
from Save_tools import filter_list 
from argparse import ArgumentParser
ROOT.gROOT.SetBatch(True)

parser = ArgumentParser()
parser.add_argument('-y', '--years', dest='years', action='store', type=str, choices=['2016', '2016APV', '2017', '2018'], default='2017')
parser.add_argument('-t', '--topology', dest='topology', action='store', type=str, default='all')
parser.add_argument('-ft', '--fittype', dest='fittype', action='store', type=str, choices=['prefit', 'postfit'], default='postfit')

args = parser.parse_args()
year = args.years
topology = args.topology
fittype = args.fittype

path_fw = os.environ['CMSSW_BASE']+"/src/MakeWorkspace_XToYYprime/"
store_path = path_fw +"/"+ year + "_" + topology

def main():

 for category in ["2fatjetsHPRest","3jetsExclRest"]:
  if "2fatjets" in category:
    mj1_bins=mj1_bins_2fatjets
    mj2_bins=mj2_bins_2fatjets
    mjj_bins=mjj_bins_2fatjets
  else:
    mj1_bins=mj1_bins_3jets
    mj2_bins=mj2_bins_3jets
    mjj_bins=mjj_bins_3jets

  colors = {
    "QCD_madgraph_pythia8_"+category: ROOT.kBlue,
    "ST_"+category: ROOT.kViolet-1,
    "TT_"+category: ROOT.kViolet,
    "WJets_"+category: ROOT.kGreen,
    "ZJets_"+category: ROOT.kPink+6,
    "VV_"+category: ROOT.kRed+2
  }
  
  names = {
      "QCD_madgraph_pythia8_"+category: "QCD",
      "ST_"+category: "Single Top",
      "TT_"+category: "t#bar{t}",
      "WJets_"+category: "W+Jets",
      "ZJets_"+category: "Z+Jets",
      "VV_"+category: "VV"
  }
  
  backgrounds = ["VV_"+category, "ST_"+category, "TT_"+category, "ZJets_"+category, "WJets_"+category, "QCD_madgraph_pythia8_"+category]
  
  samples = ["data_obs", "VV_"+category, "ST_"+category, "TT_"+category, "ZJets_"+category, "WJets_"+category, "QCD_madgraph_pythia8_"+category, "total_background"]
  
  var_list = ["fatjet","2jets","3jets"]

  hist_convert3D = {
    sample: None for sample in samples
  }
  
  hist_convert3D_proj = {
    sample : {
      var : None for var in var_list
    } for sample in samples

  }
  
  for signal_sample in signal_samples:
      signal_path = store_path + "/" + signal_sample
      os.chdir(signal_path)
      Yprime_mass =float( re.search(r'MYprime(\d+)', signal_sample).group(1) )
      print (Yprime_mass)
      Yprime_mass_up   = Yprime_mass*(1+0.3)
      Yprime_mass_down = Yprime_mass*(1-0.3)
      mj2_bins_reduce = filter_list(mj2_bins,Yprime_mass_down,Yprime_mass_up) 
      if len(mj2_bins_reduce)==0:
        mj2_bins_reduce=mj2_bins[:7]
      print (mj1_bins)
      print (mj2_bins_reduce)
      print (mjj_bins)
  
      print(signal_path+"/fitDiagnostics_CR_"+category+".root")
      f = ROOT.TFile(signal_path+"/fitDiagnostics_CR_"+category+".root")
  
      if fittype == "postfit": 
        dir_bin1 = f.Get("shapes_fit_b/c"+category)
      if fittype == "prefit": 
        dir_bin1 = f.Get("shapes_prefit/c"+category)

      n_bins_X = len(mj1_bins)-1
      n_bins_Y = len(mj2_bins_reduce)-1
      n_bins_Z = len(mjj_bins)-1
   
      data_graph = dir_bin1.Get("data")
   
      x_ctypes = ctypes.c_double()
      y_ctypes = ctypes.c_double()

      total_bkg = dir_bin1.Get("total_background")
      data_hist_1D = total_bkg.Clone("data_hist_clone_"+signal_sample)
      data_hist_1D.Reset()
   
      print("total 3D bins is ", (n_bins_X*n_bins_Y*n_bins_Z))
      print("number of data bins is ",data_hist_1D.GetNbinsX())
      

      n_points = data_graph.GetN()
      for i in range(n_points):
          data_graph.GetPoint(i, x_ctypes, y_ctypes)
          bin_num = data_hist_1D.FindBin(x_ctypes.value)
          if bin_num <= data_hist_1D.GetNbinsX() and bin_num >= 1:
              data_hist_1D.SetBinContent(bin_num, y_ctypes.value)
              data_hist_1D.SetBinError(bin_num, data_graph.GetErrorY(i))
      
      Convert1D_to_3D(mj1_bins,mj2_bins_reduce,mjj_bins, hist_convert3D, data_hist_1D, "data_obs")
     
      print(hist_convert3D)
      hist_convert3D_proj["data_obs"]["fatjet"] = hist_convert3D["data_obs"].ProjectionX("data_obs_projX")
      hist_convert3D_proj["data_obs"]["2jets"] = hist_convert3D["data_obs"].ProjectionY("data_obs_projY")
      hist_convert3D_proj["data_obs"]["3jets"] = hist_convert3D["data_obs"].ProjectionZ("data_obs_projZ")
      
      #stack = ROOT.THStack("stack", "")
    
      hists_mj1 = []
      hists_mj2 = []
      hists_mjj = []
      
      for bg in backgrounds:
          hist = dir_bin1.Get(bg)
          Convert1D_to_3D(mj1_bins,mj2_bins_reduce,mjj_bins, hist_convert3D, hist, bg)
          
          hist_convert3D[bg].SetLineColor(ROOT.kBlack)
          hist_convert3D[bg].SetFillColor(colors[bg])
          hist_convert3D[bg].SetLineWidth(1)
          hist_convert3D_proj[bg]["fatjet"] = hist_convert3D[bg].ProjectionX(bg+"_projX")
          hist_convert3D_proj[bg]["2jets"] = hist_convert3D[bg].ProjectionY(bg+"_projY")
          hist_convert3D_proj[bg]["3jets"] = hist_convert3D[bg].ProjectionZ(bg+"_projZ")
          #stack.Add(hist)
          #hists_mj1.append(hist_convert3D_proj[bg]["mj1"])
          #hists_mj2.append(hist_convert3D_proj[bg]["mj2"])
          #hists_mjj.append(hist_convert3D_proj[bg]["mjj"])
      
      total_background = dir_bin1.Get("total_background")
      Convert1D_to_3D(mj1_bins,mj2_bins_reduce,mjj_bins, hist_convert3D,total_background,"total_background")
     
      hist_convert3D_proj["total_background"]["fatjet"] = hist_convert3D["total_background"].ProjectionX("total_background_projX")
      hist_convert3D_proj["total_background"]["2jets"] = hist_convert3D["total_background"].ProjectionY("total_background_projY")
      hist_convert3D_proj["total_background"]["3jets"] = hist_convert3D["total_background"].ProjectionZ("total_background_projZ")
     
      print(hist_convert3D)
      #total_signal = dir_bin1.Get("total_signal")
      plot_data_vs_MC(fittype,year,hist_convert3D_proj,backgrounds,names,colors,"fatjet",category)
      plot_data_vs_MC(fittype,year,hist_convert3D_proj,backgrounds,names,colors,"2jets",category)
      plot_data_vs_MC(fittype,year,hist_convert3D_proj,backgrounds,names,colors,"3jets",category)

    
def plot_data_vs_MC(fittype,year,hist_convert3D_proj,backgrounds,names,colors,var,category):
 stack = ROOT.THStack("stack", "")
 for i, bkg in enumerate(backgrounds):
     print (bkg)
     hist_convert3D_proj[bkg][var].SetLineColor(ROOT.kBlack)
     hist_convert3D_proj[bkg][var].SetFillColor(colors[bkg])
     hist_convert3D_proj[bkg][var].SetLineWidth(1)
     stack.Add(hist_convert3D_proj[bkg][var])

 total_background = hist_convert3D_proj["total_background"][var]
 data_hist = hist_convert3D_proj["data_obs"][var]

 RATIO2 = total_background.Clone("RATIO_"+var)
 ratio_hist = data_hist.Clone("ratio_"+var)
 N_bins = RATIO2.GetNbinsX()
 for i in range(1, N_bins+1):
   total_bkg_content = total_background.GetBinContent(i)
   total_bkg_error = total_background.GetBinError(i)
   #print(f"bin {i} content is {total_bkg_content}")
   #print(f"bin {i} error is {total_bkg_error}")
   RATIO2.SetBinContent(i,1)
   if (total_background.GetBinContent(i)!=0):
     RATIO2.SetBinError(i,total_background.GetBinError(i)/total_background.GetBinContent(i))
     ratio_hist.SetBinContent(i,data_hist.GetBinContent(i)/total_background.GetBinContent(i))
     ratio_hist.SetBinError(i,data_hist.GetBinError(i)/total_background.GetBinContent(i))
   

 total_background.SetLineColor(ROOT.kBlack)
 total_background.SetLineWidth(2)
 total_background.SetFillStyle(0)
 
 #total_signal.SetLineColor(ROOT.kRed)
 #total_signal.SetLineWidth(3)
 #total_signal.SetLineStyle(ROOT.kDashed)
 
 data_hist.SetMarkerStyle(20)
 data_hist.SetMarkerSize(1.3)
 data_hist.SetLineColor(1)
 data_hist.SetLineWidth(2)
 data_hist.SetMarkerColor(1)
 
 canvas = ROOT.TCanvas("canvas", "Post-Fit Comparison", 800, 600)
 
 pad1 = ROOT.TPad("pad1", "pad1", 0.01, 0.32, 0.99, 0.99)
 pad1.SetBottomMargin(0.02)
 pad1.SetLeftMargin(0.11)
 pad1.SetRightMargin(0.035)
 pad1.SetTopMargin(0.08)
 pad1.Draw()
 
 pad2 = ROOT.TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.32)
 pad2.SetTopMargin(0.045)
 pad2.SetBottomMargin(0.3)
 pad2.SetLeftMargin(0.11)
 pad2.SetRightMargin(0.035)
 pad2.Draw()
 
 pad1.cd()
 pad1.SetLogy(True)  # 设置为对数坐标，如果需要线性坐标可以注释掉
 
 stack.Draw("histo")
 stack.GetYaxis().SetTitleSize(0.070)
 stack.GetXaxis().SetTitleSize(0.070)
 stack.GetYaxis().SetLabelSize(0.070)
 stack.GetXaxis().SetLabelSize(0.0)
 stack.SetTitle("")
 stack.GetYaxis().SetTitle("Events")
 stack.GetXaxis().SetTitle("")
 stack.GetYaxis().SetTitleOffset(0.80)
 stack.GetXaxis().SetTitleOffset(0.85)
 
 max_val = max(data_hist.GetMaximum(), total_background.GetMaximum())
 min_val = min([data_hist.GetBinContent(i) for i in range(1, data_hist.GetNbinsX()+1) if data_hist.GetBinContent(i) > 0] + 
               [total_background.GetBinContent(i) for i in range(1, total_background.GetNbinsX()+1) if total_background.GetBinContent(i) > 0])
 stack.SetMaximum(max_val * 100)
 #stack.SetMinimum(max(0.1, min_val * 0.1))
 stack.SetMinimum(0.5)
 
 total_background.SetFillStyle(3005)
 total_background.SetFillColor(12)  
 total_background.SetLineColor(12)  
 total_background.Draw("E2same")
 
 #total_signal.Draw("HIST SAME")
 data_hist.Draw("E same")
 
 legend = ROOT.TLegend(0.60, 0.57, 0.89, 0.91)
 legend.SetBorderSize(0)
 legend.SetFillStyle(0)
 legend.SetTextSize(0.05)
 legend.SetTextFont(62)
 
 legend.AddEntry(data_hist, "Data", "ep")
 for i, bg in enumerate(backgrounds[::-1]):
     legend.AddEntry(hist_convert3D_proj[bg][var], names[bg], "f")
 #legend.AddEntry(total_background, "Total Bkg.", "l")
 #legend.AddEntry(total_signal, "Signal", "l")
 legend.Draw()

 #cms_text = ROOT.TLatex()
 #cms_text.SetNDC()
 #cms_text.SetTextSize(0.045)
 #cms_text.DrawLatex(0.15, 0.92, "#font[62]{CMS} #font[52]{Preliminary}")
 #cms_text.DrawLatex(0.7, 0.92, "Post-Fit")
 
 pad2.cd()
 
 #ratio_hist = data_hist.Clone("ratio")
 #ratio_hist.Divide(total_background)
 
 ratio_hist.SetTitle("")
 ratio_hist.SetStats(0)
 ratio_hist.SetMarkerStyle(21)
 ratio_hist.SetMarkerSize(1.0)
 ratio_hist.SetLineColor(1)
 ratio_hist.SetMarkerColor(1)
 
 ratio_hist.GetYaxis().SetTitle("Data/Bkg.")
 ratio_hist.GetYaxis().SetNdivisions(505)
 ratio_hist.GetYaxis().SetTitleSize(0.09)
 ratio_hist.GetYaxis().SetTitleOffset(0.5)
 ratio_hist.GetYaxis().SetLabelSize(0.13)
 
 ratio_hist.GetXaxis().SetTitle(var+" mass [GeV]")
 ratio_hist.GetXaxis().SetTitleSize(0.15)
 ratio_hist.GetXaxis().SetTitleOffset(1.0)
 ratio_hist.GetXaxis().SetLabelSize(0.13)
 
 ratio_hist.SetMinimum(0.7)
 ratio_hist.SetMaximum(1.3)
 
 ratio_hist.Draw("E0psame")
 
 RATIO2.SetFillStyle(3002);
 RATIO2.SetFillColor(12);
 RATIO2.SetLineColor(12);
 RATIO2.SetMarkerSize(0);
 RATIO2.Draw("E2same");
 
 
 line = ROOT.TLine(ratio_hist.GetXaxis().GetXmin(), 1, 
                   ratio_hist.GetXaxis().GetXmax(), 1)
 line.SetLineColor(1)
 line.SetLineWidth(1)
 #line.SetLineStyle(ROOT.kDashed)
 line.Draw()
 
 
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
 #lumi_13TeV = "59.56 fb^{-1}"
 #lumiText += lumi_13TeV
 lumiText =year+" (13 TeV)"
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
 #latex.DrawLatex(l+0.12, 0.95, extraText)  
 
 
 canvas.Update()
 canvas.Draw()
 
 canvas.SaveAs(f"{fittype}_"+var+"_"+category+".pdf")

 del pad
 del pad1 
 del pad2 
 print(f"Plot saved as {fittype}_"+var+".pdf")
 
 #input("Press Enter to exit...")


def Convert1D_to_3D(X_bins,
                    Y_bins,
                    Z_bins,hist_1Dconvert3D,hist_1D,sample_name):

    print("bins", len(X_bins)-1 , len(Y_bins)-1 , len(Z_bins)-1, hist_1D.GetNbinsX())

    hist_1Dconvert3D[sample_name] = ROOT.TH3F("h3_restored_"+sample_name, "h3_restored_"+sample_name,
                        len(X_bins)-1, array('d', X_bins),
                        len(Y_bins)-1, array('d', Y_bins),
                        len(Z_bins)-1, array('d', Z_bins)
                        )

    # 遍历所有可能的(i,j,k)组合
    for i in range(1, len(X_bins)):
        for j in range(1, len(Y_bins)):
            for k in range(1, len(Z_bins)):
                # 计算对应的一维索引（与原始转换公式相同）
                index_1d = k + (len(Z_bins)-1) * ((j-1) + (len(Y_bins)-1) * (i-1))
                content = hist_1D.GetBinContent(index_1d)
                #print(i,j,k,index_1d,hist_1D.GetBinContent(index_1d))
                error = hist_1D.GetBinError(index_1d)
                #if sample_name == "total_background":
                   #print (f"totalbackground hist {index_1d} bin's content is {content}" )
                   #print (f"totalbackground hist {index_1d} bin's error is {error}" )
                hist_1Dconvert3D[sample_name].SetBinContent(i, j, k, content)
                hist_1Dconvert3D[sample_name].SetBinError(i, j, k, error)

    #hist_1Dconvert3D[sample_name] = hist_3D 
    '''
    hist_3D_mjj = hist_1Dconvert3D[sample_name].ProjectionZ("projZ") 
    hist_3D_mjj.SetTitle(sample_name+"projZ")
    c1 = ROOT.TCanvas("c1", "Combined Projections", 1200, 600)
    c1.SetLogy()
    hist_3D_mjj.Draw()
    c1.Draw()
    c1.SaveAs("convert_1Dto3D_proj/"+sample_name+".pdf")
    #del hist_3D
    del hist_3D_mjj
    del c1
    '''
if __name__ == "__main__":
  main()
