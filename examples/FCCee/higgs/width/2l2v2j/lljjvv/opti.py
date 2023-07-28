#!/usr/bin/env python

import ROOT

inputDir = 'outputs/fccee/higgs/mH-recoil/neutrinos/lljjvv/final/cutshisto/'
outputDir = 'outputs/fccee/higgs/mH-recoil/neutrinos/lljjvv/plots/opti/'

signals = ['wzp6_ee_mumuH_HZZ_ecm240',
           'wzp6_ee_eeH_HZZ_ecm240']
backgrounds = ['wzp6_ee_mumuH_HWW_ecm240',
               'wzp6_ee_mumuH_HZa_ecm240',
               'wzp6_ee_mumuH_Haa_ecm240',
               'wzp6_ee_mumuH_Hbb_ecm240',
               'wzp6_ee_mumuH_Hcc_ecm240',
               'wzp6_ee_mumuH_Hgg_ecm240',
               'wzp6_ee_mumuH_Hmumu_ecm240',
               'wzp6_ee_mumuH_Hss_ecm240',
               'wzp6_ee_mumuH_Htautau_ecm240',
               'wzp6_ee_nunuH_HZZ_ecm240',
               'wzp6_ee_nunuH_HWW_ecm240',
               'wzp6_ee_nunuH_HZa_ecm240',
               'wzp6_ee_nunuH_Hbb_ecm240',
               'wzp6_ee_nunuH_Hcc_ecm240',
               'wzp6_ee_nunuH_Hgg_ecm240',
               'wzp6_ee_nunuH_Hmumu_ecm240',
               'wzp6_ee_nunuH_Hss_ecm240',
               'wzp6_ee_nunuH_Htautau_ecm240',
               #'wzp6_ee_mumu_ecm240',
               'wzp6_ee_eeH_HWW_ecm240',
               'wzp6_ee_eeH_HZa_ecm240',
               'wzp6_ee_eeH_Haa_ecm240',
               'wzp6_ee_eeH_Hbb_ecm240',
               'wzp6_ee_eeH_Hcc_ecm240',
               'wzp6_ee_eeH_Hgg_ecm240',
               'wzp6_ee_eeH_Hmumu_ecm240',
               'wzp6_ee_eeH_Hss_ecm240',
               'wzp6_ee_eeH_Htautau_ecm240',
               #'wzp6_ee_ee_Mee_30_150_ecm240',
               'p8_ee_ZZ_ecm240',
               'p8_ee_WW_ecm240']

var = ['leptonic_recoil_m'] #par exemple, mais osef
sels = ['sel0', 
	'sel1',
	'sel2',
	'sel3', 
	'sel4', 
	'sel5', 
	'sel6', 
	'sel7', 
	'sel8', 
	'sel9', 	
	'sel10',
	'sel11' 	
#	'sel12', 	
#	'sel13',	
#	'sel14'	
	]

#cuts = []
#cut = 78
#cuts.append(f"{cut}") 
#for j in range (1,11):
#    cut = cut + 5
#    cuts.append(f"{cut}")
#
#print(cuts)  

cuts = []
cut = 128
cuts.append(f"{cut}") 
for j in range (1,11):
    cut = cut + 1
    cuts.append(f"{cut}")

print(cuts)              


formats = ['png']
final = 0 
def main():
    n = len(sels)
    print(n)
    h = ROOT.TH1F("cutssoversqrtb", "S over sqrtB as a function of the value of the cut on recoil mass", n+1, 0, n+1)
    xax = h.GetXaxis()
    for i in range(1,n):
        print(i)
        sums = 0 
        sumb = 0 
        for s in signals:
                fs = ROOT.TFile.Open(f"{inputDir}{s}_sel{i}_histo.root")
                for vs in var:
                    hs = fs.Get(vs) 
                sums = sums + hs.Integral()
                #print(sums*5*10**(6))
                
        for b in backgrounds: 
                fb = ROOT.TFile.Open(f"{inputDir}{b}_sel{i}_histo.root")
                for vb in var: 
                    hb = fb.Get(vb) 
                sumb = sumb + hb.Integral()
                
      
        final = (sums/(sumb**(1/2)))*(5*10**(6))**(1/2)
        h.SetBinContent(i, final) 
        xax.SetBinLabel(i, cuts[i])
        

    c = ROOT.TCanvas()
    h.SetXTitle(f"Value of the cut on {var[0]}") 
        
    #xax = ROOT.TGaxis(15,0,65,3,"h", 10)
    #xax.Draw()
    
    #xax = h.GetXaxis()
    #xax.SetLimits(0,18)	
    h.Draw("hist")
    #axis8 = ROOT.TGaxis(8,-0.8,8,0.8,0,9000,50510,"+L");
    #axis8.Draw()
    for ext in formats:
        c.SaveAs(f"{outputDir}{var[0]}.{ext}")
    
     

    print(final) 
       

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    main()

