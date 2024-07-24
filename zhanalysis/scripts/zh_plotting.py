import ROOT 
import numpy as np
import math 
#from zhanalysis.stage1.ZH_Hadronic_stage1 import vars


algs = ["anti-kt", "durham-kt"]
#select algorithm for plot 0--anti-kt  1--durham-kt
alg = 0

#read ROOT file
#file = ROOT.TFile("chunk_1.root", "READ")
file = ROOT.TFile("ccH_Hbb_1620_10k.root", "READ")

#file = ROOT.TFile("chunk_0.root", "READ")
#get tree in Root File
antikt_tree = file.Get("events")
print("finished reading")


#vector with c jet pair invariant masses
c_pairs = []
#vector with b jet pair invariant masses
b_pairs = []

counter = 0

for event in antikt_tree:

   if counter >= 5000:
      counter+=1
      break  # Stop after processing 5000 events
      

   if (event.event_njet==4):
         counter+=1
         print(counter)
         print(event.jets_truth)
         #events4.append(event.event_njet)
         c_idx=np.where(np.abs(event.jets_truth)==4)[0]
         b_idx=np.where(np.abs(event.jets_truth)==5)[0]
         

         if (len(c_idx)==2 and len(b_idx)==2):
         
            idx_map = {
            "c": c_idx,
            "b": b_idx}

            pairs_map = {
            "c": c_pairs,
            "b": b_pairs    
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
               # vec1.SetPxPyPzE(event.jet_px_corr[int(idx[0])],event.jet_py_corr[int(idx[0])],event.jet_pz_corr[int(idx[0])],event.jet_e_corr[int(idx[0])])
               # vec2.SetPxPyPzE(event.jet_px_corr[int(idx[1])],event.jet_py_corr[int(idx[1])],event.jet_pz_corr[int(idx[1])],event.jet_e_corr[int(idx[1])])
               # #initialize sum of pair tlv 
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

#if counter>=4000:
   #compute number of bins and xmax and xmin for c masses and b masses
   max_cc = int(max(c_pairs))
   min_cc = int((min(c_pairs)))
   bins_cc = int(math.sqrt(len(c_pairs)))

   max_bb = int(max(b_pairs))
   min_bb = int(min(b_pairs))
   bins_bb = int(math.sqrt(len(b_pairs)))

   # print(max_cc)
   # print(max_bb)

   #define number of bins and xmax and xmin for bb jets and cc jets
   nbins_cc = bins_cc
   xmin_cc = min_cc
   xmax_cc = max_cc

   nbins_bb = bins_bb
   xmin_bb = min_bb
   xmax_bb = max_bb

   print("done reading")
   "hist_m", "Jet Mass; Jet Mass; Entries"
   # #create histogram from TH1 class
   hist_cc = ROOT.TH1F("cc Jet Pair Mass", "{} Algorithm".format(algs[alg]), nbins_cc, xmin_cc, xmax_cc)
   hist_bb =  ROOT.TH1F("bb Jet Pair Mass", "{} Algorithm".format(algs[alg]), nbins_bb, xmin_bb, xmax_bb)

   print("done initializing")

   for i in range(len(c_pairs)):
      cc = c_pairs[i]
      bb = b_pairs[i]
      hist_cc.Fill(cc)
      hist_bb.Fill(bb)

   print("histograms filled")

   #create canvas class --- should produce plots
   canvas1 = ROOT.TCanvas("canvas1", "cc and bb Jet masses", 800, 600)

   hist_cc.SetLineColor(ROOT.kBlue)
   hist_bb.SetLineColor(ROOT.kRed)

   hist_cc.Draw()
   hist_bb.Draw("SAME")

   legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
   legend.AddEntry(hist_cc, "cc Jet Mass", "l")
   legend.AddEntry(hist_bb, "bb Jet Mass", "l")
   legend.Draw()

   canvas1.SaveAs("durhamtest/durhamccbb.pdf")


   canvas2 = ROOT.TCanvas("canvas2", "bb Jet masses", 800, 600)

