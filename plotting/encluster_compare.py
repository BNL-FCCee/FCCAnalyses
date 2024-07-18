import ROOT 
import numpy as np
import math 

#returns a vector with cpair vector as first entry and bpair vector as second entry
colors=[ROOT.kMagenta, ROOT.kRed, ROOT.kBlue, ROOT.kOrange, ROOT.kGreen, ROOT.kCyan]

reco_map={}

#call files
#files = ["ccH_Hbb_1620_10k.root", "ccH_Hbb_410_10k.root","chunk_1.root"]
files = ["ccH_Hbb_1620_10k.root","chunk_1.root"]
#set number of events
nevents = 4000
#set number of bins
bins = 35

## 0--indicates recojet_px etc. , 1-- indicates jet_px_corr etc.
#algos = [0,0,1]
algos = [0,1]
#write labels associated with file in order 
algs = ["anti-kt r=1.6 e-cut=20GeV", "anti-kt r=0.4 e-cut=10GeV","durham-kt with correction","cambridge"]
algs = ["anti-kt r=1.6 e-cut=20GeV", "durham-kt with correction","cambridge"]


for i in range(len(files)):
    reco_map[files[i]]=algos[i]

#initialize pair lists
c_pairs =[]
b_pairs = []


def get_flav_idx(p_idx, event):
    jets_e = np.array(event.jet_e)
    #choose jet energies indexed by p_idx
    p_e = jets_e[p_idx]

    #create array with energies in order of highest to lowest
    p_hilo = sorted(p_e, reverse=True)
    
    #create array with jet energy indices in order of highest to lowest energy
    p_e_idx = []
    for p in p_hilo:
        id = np.where(p_e==p)[0][0]
        p_e_idx.append(id)

    #select the two highest indices to represent the b and c candidates
    p_idx_2 = p_e_idx[:2]

    return p_idx_2

def get_indices(event,flav)-> list:
    truth_list = np.array(event.jets_truth)
    c_idx=np.where(np.abs(truth_list)==4)[0]
    b_idx=np.where(np.abs(truth_list)==5)[0]
    clen = len(c_idx) 
    # print("clen",clen)
    blen = len(b_idx)
    # print("blen",blen)


    if (clen>=2 and blen>=2): 
        idx=[]
        if not (clen==2 and blen==2):
            #find indices of jet energies in order of energy
            c_idx = get_flav_idx(c_idx,event)
            # print("cidx after",len(c_idx))
            b_idx = get_flav_idx(b_idx,event)
            # print("bidx after",len(b_idx))
        if flav==0:
            for c in c_idx:
                idx.append(c)
        if flav==1:
            for b in b_idx:
                idx.append(b)
        return idx

def create_vectors(alg,idx,event):
    #set tlv values -- 
    vec1 = ROOT.TLorentzVector() 
    vec2 = ROOT.TLorentzVector()
    vecs = [vec1,vec2]
    for i, vec in enumerate(vecs):

        if alg==0:
            #set tlv values-- Anti-kt Alg
            vec.SetPxPyPzE(event.recojet_px[int(idx[i])],event.recojet_py[int(idx[i])],event.recojet_pz[int(idx[i])],event.jet_e[int(idx[i])])
            
        if alg ==1:
            #set tlv values-- Durham Alg
            vec.SetPxPyPzE(event.jet_px_corr[int(idx[i])],event.jet_py_corr[int(idx[i])],event.jet_pz_corr[int(idx[i])],event.jet_e_corr[int(idx[i])])

    return vecs

#flav-- 0 is c  flav-- 1 is b
def create_pair(file, flav):

    pairs = []

    #determine which reco jet variable to use
    alg = reco_map[file]

    file =  ROOT.TFile(file, "READ")
    tree = file.Get("events")

    print("finished reading")

    for i, event in enumerate(tree):
        if i >= nevents:  # Stop after nevents
            break
        if all(truth is not None for truth in event.jets_truth):
            if event.event_njet>=4:
                # print(i)

                idx=get_indices(event,flav)

                if idx is not None:
                
                    #create two vectors
                    vectors = create_vectors(alg,idx,event)
                    vec1=vectors[0]
                    vec2=vectors[1]
                    
                    #initialize sum of pair tlv 
                    sum = ROOT.TLorentzVector()

                    #compute tlv sum
                    sum = vec1 + vec2
                    
                    #get invariant mass of the sum 
                    mass=sum.M()
                    
                    #add invariant mass to pairs list
                    pairs.append(mass)
            
    return pairs

#pairs is a list of c or b pair lists 
#flav-- 0 is c  flav-- 1 is b
def create_hist(pairs, flav):

    maximum = 160
    minimum = 0

    flavs = ["cc","bb"]

    hists = []
    for i, pp in enumerate(pairs):
        
        hist = ROOT.TH1F("hist"+str(i), flavs[flav]+" jet pair masses", bins, minimum, maximum)
    
        for p in pp:
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
        hist.SetMaximum(hist.GetMaximum() * 1.1) 
        hist.Sumw2(ROOT.kFALSE)
        if index==0:
            hist.Draw("HIST")
        else:
            hist.Draw("SAME")
        print(algs[index]+"mean")
        print(hist.GetMean())
    legend.SetBorderSize(0)
    legend.Draw()
    canvas.SaveAs("antikt_4jetse/10k"+flavs[flav]+"mass.pdf")

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

