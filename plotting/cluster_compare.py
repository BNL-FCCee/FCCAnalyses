import ROOT 
import numpy as np
import math 
#from zhanalysis.stage1.ZH_Hadronic_stage1 import vars

   
#returns a vector with cpair vector as first entry and bpair vector as second entry
colors=[ROOT.kOrange, ROOT.kMagenta, ROOT.kBlue, ROOT.kRed, ROOT.kGreen]

#files = ["ccH_Hbb_1620_10k.root"]
files = ["ccH_Hbb_1620_10k.root","chunk_1.root"]

#write labels associated with file in order 
algs = ["anti-kt, r-1.6, ecut-20", "durham-kt","cambridge"]

c_pairs =[]
b_pairs = []


#flav-- 0 is c  flav-- 1 is b
def create_pair(file, flav):

    pairs = []
    
    file =  ROOT.TFile(file, "READ")
    tree = file.Get("events")
    print("finished reading")
    for i, event in enumerate(tree):
        if i >= 100:  # Stop after 10,000 events
            break
        if (event.event_njet==4):
            #print(i)
            c_idx=np.where(np.abs(event.jets_truth)==4)[0]
           # print(c_idx)
            b_idx=np.where(np.abs(event.jets_truth)==5)[0]
           
            if (len(c_idx)==2 and len(b_idx)==2):
            
                if flav==0:
                    idx = c_idx
                if flav==1:
                    idx = b_idx
                
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
                
                # print(mass)
                #add invariant mass to 
                pairs.append(mass)
    return pairs


#pairs is a vector of c or b pair vectors 
#flav-- 0 is c  flav-- 1 is b
def create_hist(pairs, flav):

    maximum = pairs[0][0]
    minimum = pairs[0][0]
    bins = int(math.sqrt(int(len(pairs[0]))))

    
    flavs = ["cc","bb"]

    hists = []
    for i, pp in enumerate(pairs):
        maxi = max(pp)
        if maxi>maximum:
            maximum = maxi
        mini = min(pp)
        if mini<minimum:
            minimum = mini
        hist = ROOT.TH1F("hist"+str(i), flavs[flav]+" jet pair masses", bins, maximum, minimum)
        hist.Scale(1.0 / hist.Integral())
        for p in pp:
            print(p)
            hist.Fill(p)
        hists.append(hist)
       

    canvas = ROOT.TCanvas("canvas", flavs[flav]+" jet pair masses", 800, 600)
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    
    for index, hist in enumerate(hists):
        color = colors[index]
        hist.SetLineColor(color)
        legend.AddEntry(hist, algs[index], "l")
        if index==0:
            hist.Draw()
        else:
            hist.Draw("SAME")

    legend.Draw()
    canvas.SaveAs("1620"+flavs[flav]+"mass.pdf")



#Now run functions

for file in files:
   #create vector of c_pair vectors for each file
   f1= create_pair(file,0)
   f2= create_pair(file,1)
   c_pairs.append(f1)
   #create vector of b_pair vectors for each file
   b_pairs.append(f2)

#create cc and bb histograms
create_hist(c_pairs, 0)
create_hist(b_pairs, 1)

