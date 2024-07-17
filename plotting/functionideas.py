import ROOT 
import numpy as np
import math 
#from zhanalysis.stage1.ZH_Hadronic_stage1 import vars

def create_pairs(self, file, cpairs, bpairs):
    self.file = file
    self.cpairs = cpairs
    self.bpairs = bpairs
    cpairs = []
    bpairs = []
    pairs = []
    file =  ROOT.TFile(file, "READ")
    tree = file.Get("events")

    for event in tree:
        if (event.event_njet==4):
            c_idx=np.where(np.abs(event.jets_truth)==4)[0]
            b_idx=np.where(np.abs(event.jets_truth)==5)[0]

            if (len(c_idx)==2 and len(b_idx)==2):
            
                idx_map = {
                "c": c_idx,
                "b": b_idx}

                pairs_map = {
                "c": cpairs,
                "b": bpairs    
                }

                flavors = ["c","b"]

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
    pairs.append(cpairs)
    pairs.append(bpairs)
    return pairs

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
    
        for p in pp:
            print(p)
            hist.Fill(p)

        print(hist.GetEntries())
        hist.Scale(1.0 / hist.Integral())
        hists.append(hist)
       
import ROOT 
import numpy as np
import math 
#from zhanalysis.stage1.ZH_Hadronic_stage1 import vars

   
#returns a vector with cpair vector as first entry and bpair vector as second entry
colors=[ROOT.kMagenta, ROOT.kRed, ROOT.kBlue, ROOT.kOrange, ROOT.kGreen, ROOT.kCyan]

#files = ["ccH_Hbb_1620_10k.root"]
file1 = "ccH_Hbb_1620_10k.root"
file2 = "chunk_1.root"
files = [file1,file2]

## 0--indicates recojet_px etc. , 1-- indicates jet_px_corr etc.
reco_map = {file1: 0,
            file2: 1}

#write labels associated with file in order 
algs = ["anti-kt r=1.6 Ecut=20GeV", "durham-kt","cambridge"]
algos = ["anti-kt", "durham-kt","cambridge"]
c_pairs =[]
b_pairs = []


def create_vectors(alg,idx,event):
    #set tlv values -- 
    #initialize tlvs for c pair and b pair
    vec1 = ROOT.TLorentzVector() 
    vec2 = ROOT.TLorentzVector()

    if alg==0:
        #set tlv values-- Anti-kt Alg
        vec1.SetPxPyPzE(event.recojet_px[int(idx[0])],event.recojet_py[int(idx[0])],event.recojet_pz[int(idx[0])],event.jet_e[int(idx[0])])
        vec2.SetPxPyPzE(event.recojet_px[int(idx[1])],event.recojet_py[int(idx[1])],event.recojet_pz[int(idx[1])],event.jet_e[int(idx[1])])
    if alg ==1:
        #set tlv values-- Durham Alg
        vec1.SetPxPyPzE(event.jet_px_corr[int(idx[0])],event.jet_py_corr[int(idx[0])],event.jet_pz_corr[int(idx[0])],event.jet_e_corr[int(idx[0])])
        vec2.SetPxPyPzE(event.jet_px_corr[int(idx[1])],event.jet_py_corr[int(idx[1])],event.jet_pz_corr[int(idx[1])],event.jet_e_corr[int(idx[1])])

    vecs = [vec1, vec2]
    return vecs



#flav-- 0 is c  flav-- 1 is b
def create_pair(file, flav):

    pairs = []
    
    file =  ROOT.TFile(file, "READ")
    tree = file.Get("events")
    print("finished reading")
    for i, event in enumerate(tree):
        if i >= 2000:  # Stop after 10,000 events
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

                alg = reco_map[file]

                vecs = create_vectors(alg,idx,event)
                vec1=vecs[0]
                vec2=vecs[1]
                

                # #initialize tlvs for c pair and b pair
                # vec1 = ROOT.TLorentzVector() 
                # vec2 = ROOT.TLorentzVector()

                # #set tlv values -- 
                # vec1.SetPxPyPzE(event.recojet_px[int(idx[0])],event.recojet_py[int(idx[0])],event.recojet_pz[int(idx[0])],event.jet_e[int(idx[0])])
                # vec2.SetPxPyPzE(event.recojet_px[int(idx[1])],event.recojet_py[int(idx[1])],event.recojet_pz[int(idx[1])],event.jet_e[int(idx[1])])
                
                # #set tlv values-- Durham Alg
                # vec1.SetPxPyPzE(event.recojet_px[int(idx[0])],event.recojet_py[int(idx[0])],event.recojet_pz[int(idx[0])],event.jet_e[int(idx[0])])
                # vec2.SetPxPyPzE(event.recojet_px[int(idx[1])],event.recojet_py[int(idx[1])],event.recojet_pz[int(idx[1])],event.jet_e[int(idx[1])])


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

    maximum = 160
    minimum = 0
    bins = 32

    flavs = ["cc","bb"]

    hists = []
    for i, pp in enumerate(pairs):
        
        hist = ROOT.TH1F("hist"+str(i), flavs[flav]+" jet pair masses", bins, minimum, maximum)
    
        for p in pp:
            print(p)
            hist.Fill(p)

        print(hist.GetEntries())
        hist.Scale(1.0 / hist.Integral())
        hists.append(hist)
       

    canvas = ROOT.TCanvas("canvas", flavs[flav]+" jet pair masses", 800, 600)
    legend = ROOT.TLegend(0.4, 0.5, 0.1, 0.7)
    
    for index, hist in enumerate(hists):
        color = colors[index]
        hist.SetLineColor(color)
        legend.AddEntry(hist, algs[index], "l")
        hist.Sumw2(ROOT.kFALSE)
        if index==0:
            hist.Draw("HIST")
        else:
            hist.Draw("SAME")
        print(algs[index])
        print(hist.GetMean())
    legend.SetBorderSize(0)
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

