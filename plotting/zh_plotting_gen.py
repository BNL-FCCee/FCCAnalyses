import ROOT 
import numpy as np
import math 
#from zhanalysis.stage1.ZH_Hadronic_stage1 import vars

in_files = ["ccH_Hbb_1620_10k.root","ccH_Hbb_410k.root"]
filelen = len(in_files)

#vector of root files
files = []
#vector of trees
trees=[]

#vector with c jet pair invariant masses
c_pairs = []
#vector with b jet pair invariant masses
b_pairs = []


for i in range(filelen):
#declare ROOT files
    files.append(ROOT.TFile(in_files[i], "READ"))
    trees.append(files[i].Get("events"))
    c_pairs.append([])  
    b_pairs.append([])

for event in trees[i]:
   #create vector that will contain events with 4 jets
   if (event.event_njet==4):

      #events4.append(event.event_njet)
      c_idx=np.where(np.abs(event.jets_truth)==4)[0]
      b_idx=np.where(np.abs(event.jets_truth)==5)[0]

      if (len(c_idx)==2 and len(b_idx)==2):
      
         idx_map = {
         "c": c_idx,
         "b": b_idx}

         pairs_map = {
         "c": c_pairs[i],
         "b": b_pairs[i]    
         }

         flavors = ["c","b"]
         vec_d={}

         for flav in flavors:
            #define index variable from dict.
            idx = idx_map[flav]
            #define pairs variable from dict. 
            pairs = pairs_map[flav]

            #initialize tlvs for c pair and b pair
            vec1 = ROOT.TLorentzVector() 
            vec2 = ROOT.TLorentzVector()

            #set tlv values
            vec1.SetPxPyPzE(event.recojet_px[int(idx[0])],event.recojet_py[int(idx[0])],event.recojet_pz[int(idx[0])],event.jet_e[int(idx[0])])
            vec2.SetPxPyPzE(event.recojet_px[int(idx[1])],event.recojet_py[int(idx[1])],event.recojet_pz[int(idx[1])],event.jet_e[int(idx[1])])
            
            #initialize sum of pair tlv 
            sum = ROOT.TLorentzVector()

            #compute tlv sum
            sum = vec1 + vec2

            #get invariant mass of the sum 
            mass=sum.M()

            #to refer to invariant mass of c pair and b pair
            vec_d[flav]= mass
            # print(mass)
            #add invariant mass to 
            pairs.append(mass)

#define dictionaries for max, min, number of bins for cc and bb jets


hist_ccd={}
hist_bbd={}

for c_pair in c_pairs:
    
    nbins_cc = int(max(c_pair))
    xmin_cc  = int(min(c_pair))
    xmax_cc = int(math.sqrt(len(c_pair)))
    hist_ccd[c_pair] = ROOT.TH1F("B Jet Pair Mass", "Anti-kt alg", nbins_cc, xmin_cc, xmax_cc)

    for cc in c_pair:
        hist_ccd[c_pair].Fill(cc)

for b_pair in b_pairs:

    nbins_bb = int(max(b_pair))
    xmin_bb  = int(min(b_pair))
    xmax_bb = int(math.sqrt(len(c_pair)))
    hist_bbd[b_pair] = ROOT.TH1F(" Jet Pair Mass", "Anti-kt alg", nbins_cc, xmin_cc, xmax_cc)

    for bb in b_pair:
        hist_bbd[c_pair].Fill(bb)

print("histograms filled")


#create canvas class --- should produce plots
canvas1 = ROOT.TCanvas("canvas1", "cc  and bb Jet masses", 800, 600)


for c_pair in c_pairs:
    
    hist_ccd[c_pair].SetLineColor(ROOT.kBlue)

    hist_cc.Draw()
    hist_bb.Draw("SAME")

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(hist_cc, "cc Jet Mass", "l")
legend.AddEntry(hist_bb, "bb Jet Mass", "l")
legend.Draw()

canvas1.SaveAs("antikt1620ccbbmass.pdf")


canvas2 = ROOT.TCanvas("canvas2", "bb Jet masses", 800, 600)

       
      