import sys
import os
import glob
import string
import subprocess
import ROOT
import re
from array import array
#from commands import getoutput
from parameter import signal_samples
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-y', '--years', dest='years', action='store', type=str, choices=['2016', '2016APV', '2017', '2018'], default='2017')
parser.add_argument('-t', '--topology', dest='topology', action='store', type=str, choices=['boosted', 'resolved'], default='resolved')


args = parser.parse_args()
year = args.years
topology = args.topology

#mass = [re.search(r'MX\d+', s).group() for s in signal_samples]
Xmass = [int(re.search(r'MX(\d+)', s).group(1)) for s in signal_samples] 

xsec=[1, 1, 1]

limit_Obs=[]

path_fw = os.environ['CMSSW_BASE']+"/src/MakeWorkspace_XToYYprime/"
store_path = path_fw +"/"+ year + "_" + topology

def main():
  os.chdir(store_path)
  outputpath = path_fw 
  limit_Exp_p2s=[]
  limit_Exp_p1s=[]
  limit_Exp=[]
  limit_Exp_m1s=[]
  limit_Exp_m2s=[]
  mass_err=[]
  Xmass_rm_err=Xmass[:]
  xsec_rm_err=xsec[:]

  readCombine(store_path,signal_samples,limit_Exp_m2s, limit_Exp_m1s, limit_Exp, limit_Exp_p1s, limit_Exp_p2s,mass_err)
  
  for err in mass_err:
    index = mass.index(err)
    Xmass_rm_err.remove(Xmass[index])
    xsec_rm_err.remove(xsec[index])
      
  #if len(Xmass_rm_err) < 1:
  #  pass

  limit_Exp_m2s= [x * y for x, y in zip(xsec_rm_err, limit_Exp_m2s)]
  limit_Exp_m1s= [x * y for x, y in zip(xsec_rm_err, limit_Exp_m1s)]
  limit_Exp= [x * y for x, y in zip(xsec_rm_err, limit_Exp)]
  limit_Exp_p1s= [x * y for x, y in zip(xsec_rm_err, limit_Exp_p1s)]
  limit_Exp_p2s= [x * y for x, y in zip(xsec_rm_err, limit_Exp_p2s)]

  plot_limit(Xmass_rm_err,xsec_rm_err, limit_Exp_m2s, limit_Exp_m1s, limit_Exp, limit_Exp_p1s, limit_Exp_p2s)

def readCombine(store_path,signal_samples,Limit_Exp_m2s, Limit_Exp_m1s, Limit_Exp, Limit_Exp_p1s, Limit_Exp_p2s,Mass_err):
  for signal_sample in signal_samples:
    xmass = re.search(r'MX\d+', signal_sample).group()
    Input = store_path+"/"+signal_sample + "/combine_result_"+xmass+ ".txt"
    limitcal = False
    with open(Input,"r") as file:
      for line_number, line in enumerate(file,1):
        if "Expected  2.5%" in line:
            r_m2s = float(line.split()[4])
            #print r_m2s
            Limit_Exp_m2s.append(r_m2s)
            limitcal = True
        if "Expected 16.0%" in line:
            r_m1s = float(line.split()[4])
            #print r_m1s
            Limit_Exp_m1s.append(r_m1s)
            limitcal = True
        if "Expected 50.0%" in line:
            r = float(line.split()[4])
            #print r
            Limit_Exp.append(r)
        if "Expected 84.0%" in line:
            r_p1s = float(line.split()[4])
            #print r_p1s
            Limit_Exp_p1s.append(r_p1s)
            limitcal = True
        if "Expected 97.5%" in line:
            r_p2s = float(line.split()[4])
            #print r_p2s
            Limit_Exp_p2s.append(r_p2s)
            limitcal = True

      if not limitcal:
        Mass_err.append(Masspoint)
        print ("the mass point %s combine calculation has error" % (Masspoint))

    print ("Limit_Exp_m2s=", Limit_Exp_m2s)
    print ("Limit_Exp_m1s=", Limit_Exp_m1s)
    print ("Limit_Exp=", Limit_Exp)
    print ("Limit_Exp_p1s=", Limit_Exp_p1s)
    print ("Limit_Exp_p2s=", Limit_Exp_p2s)

def plot_limit(XMass,Xsec, Limit_Exp_m2s, Limit_Exp_m1s, Limit_Exp, Limit_Exp_p1s, Limit_Exp_p2s):
    XMass_reverse = XMass[::-1]
    Limit_Exp_p1s_reverse = Limit_Exp_p1s[::-1]
    Limit_Exp_p2s_reverse = Limit_Exp_p2s[::-1]

    XMass_sigma = XMass + XMass_reverse
    XMass_sigma.append(XMass[0])
    Limit_1sigma = Limit_Exp_m1s + Limit_Exp_p1s_reverse    
    Limit_1sigma.append(Limit_Exp_m1s[0])
    Limit_2sigma = Limit_Exp_m2s + Limit_Exp_p2s_reverse     
    Limit_2sigma.append(Limit_Exp_m2s[0])

    canvas = ROOT.TCanvas("canvas", "Combine Limits", 800, 600)
    canvas.SetLogy()
    limit_exp = ROOT.TGraph(len(XMass), array('d', XMass), array('d', Limit_Exp))
    XSEC = ROOT.TGraph(len(XMass), array('d', XMass), array('d', Xsec))
    band_1s = ROOT.TGraph(len(XMass_sigma), array('d', XMass_sigma), array('d', Limit_1sigma))
    band_2s = ROOT.TGraph(len(XMass_sigma), array('d', XMass_sigma), array('d', Limit_2sigma))

    band_2s.SetLineStyle(1);
    band_2s.SetLineColor(0);
    band_2s.SetFillColor(ROOT.kOrange);
    band_2s.Draw("AF");
    band_2s.SetTitle("");
    band_2s.GetXaxis().SetTitle("#it{m}_{X} [GeV]");
    band_2s.GetYaxis().SetTitle("#sigma(pp #rightarrow X)[fb]");
    band_2s.SetMaximum(100);

    band_1s.SetLineStyle(1);
    band_1s.SetLineColor(0);
    band_1s.SetFillColor(ROOT.kGreen+1);
    band_1s.Draw("LF2");

    limit_exp.SetLineWidth(2);
    limit_exp.SetLineStyle(2);
    limit_exp.SetMarkerSize(1.3);
    limit_exp.Draw("L");

    XSEC.SetLineColor(4);
    XSEC.SetMarkerColor(4);
    XSEC.SetLineWidth(2);
    XSEC.SetLineStyle(9);
    #XSEC.Draw("L");

    pad = ROOT.TPad("pad","pad",0.01,0.01,0.99,0.99)
    ROOT.gPad.RedrawAxis()
    #channelText = ""
    #widthText = "#Gamma/#it{m}_{Z'}="+WIDTH[0]+"%"
    #channelTextFont   = 42
    #channelTextSize   = 0.060
    cmsText     = "CMS"
    cmsTextFont   = 61  # default is helvetic-bold
    #writeExtraText = true
    #extraText   = "Simulation"
    #extraTextFont = 52  # default is helvetica-italics
    #extraText = widthText + ' '+ channelText
    extraTextFont = 42
    extraTextSize = 0.045
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
    lumi_13TeV = "59.83 fb^{-1}"
    #lumiText += lumi_13TeV
    lumiText = "2018 (13 TeV)"
    t = pad.GetTopMargin()
    b = pad.GetBottomMargin()
    r = pad.GetRightMargin()
    l = pad.GetLeftMargin()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    #extraTextSize = extraOverCmsTextSize*cmsTextSize
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(lumiTextSize*t)
    latex.DrawLatex(1-r,0.92,lumiText)
    latex.SetTextFont(cmsTextFont)
    latex.SetTextAlign(11)
    latex.SetTextSize(cmsTextSize*t)
    latex.DrawLatex(l, 0.92,cmsText)#0.075
    latex.SetTextFont(extraTextFont)
    latex.SetTextSize(extraTextSize)
    #latex.DrawLatex(0.51, 0.57, extraText)
    #latex.DrawLatex(l+0.18, 1-t+lumiTextOffset*t-0.12,extraText)#0.075
    del pad

    #widthlegend = "Theoretical #sigma, #Gamma/#it{m}_{Z'}="+WIDTH[0]+"%"
    legend = ROOT.TLegend(0.50, 0.63, 0.88, 0.88, "95% CL upper limits")
    legend.SetTextSize(0.038)
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.AddEntry(limit_exp, "Median expected", "L")
    legend.AddEntry(band_1s, "68% expected", "f")
    legend.AddEntry(band_2s, "95% expected", "f")
    #legend.AddEntry(XSEC,widthlegend,"L");
    legend.Draw()

    #canvas.Draw()
    canvas.SaveAs("./limitplot.pdf")

if __name__=='__main__':
 main()
