import uproot
import ROOT
import numpy as np
import awkward as ak


colors= [ROOT.kAzure+6, ROOT.kPink+9, ROOT.kGreen+2,
ROOT.kOrange+3, ROOT.kMagenta+2,ROOT.kCyan+2,]

files =  "hdecay.root"

file = uproot.open(files)
print("read file")
tree = file['events']

branches = tree.arrays()

#branches = tree.arrays(library = "np")
print("branches retrieved")

event_njet = branches["event_njet"]

jets_mask = event_njet==4
mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)
mask_c = np.abs(branches["jets_truth"][jets_mask][mask])==4
mask_b = np.abs(branches["jets_truth"][jets_mask][mask])==5 



# mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
# mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)
# mask_c = np.abs(branches["jets_truth"][jets_mask][mask])==4
# mask_b = np.abs(branches["jets_truth"][jets_mask][mask])==5 


# jc_phi = branches["reco_jc_phi"] 
# jc_theta = branches["reco_jc_theta"] 
# jc_e = branches["reco_jc_e"] 

#extract data for 
px = branches["jet_px_corr"]
py = branches["jet_py_corr"]
pz = branches["jet_pz_corr"]
e = branches["jet_e_corr"]
data = [px,py,pz,e]

c_data = []
b_data = []
for d in data:
    c_data.append(d[jets_mask][mask][mask_c])
    b_data.append(d[jets_mask][mask][mask_b])

print("c and b data sorted")
#quark combinations

# def get_deltaR(data):
#     for i in len(data[0]):
#         p1 = ROOT.TLorentzVector()
#         p2 = ROOT.TLorentzVector()
#         p1.set(data[0][i][0],data[1][i][0],data[2][i][0],data[3][i][0])
#         p2.set(data[0][i][1],data[1][i][1],data[2][i][1],data[3][i][1])
#         deltaR = p1.DeltaR(p2)

def set_particle_tlv(p_data):
    p_tlv = []
    for i in range(len(p_data[0])):  
        pair = []
        p1 = ROOT.TLorentzVector()
        p2 = ROOT.TLorentzVector() 
        p1.SetPxPyPzE(data[0][i][0],data[1][i][0],data[2][i][0],data[3][i][0])
        p2.SetPxPyPzE(data[0][i][1],data[1][i][1],data[2][i][1],data[3][i][1])  
        pair.append(p1)
        pair.append(p2)
        p_tlv.append(pair)
    print("particle tlvs are set")
    return p_tlv

def get_deltaR(l1,l2):
    deltaR = l1.DeltaR(l2)
    return deltaR


c_tlv = set_particle_tlv(c_data)
b_tlv = set_particle_tlv(b_data)

print("pairs of p1p2 lorentz vectors are in respective lists")

c1 = []
c2 = [] 
b1 = [] 
b2 = [] 
 
for i in range(len(c_tlv)):
    c1.append(c_tlv[i][0])
    c2.append(c_tlv[i][1])
    b1.append(b_tlv[i][0])
    b2.append(b_tlv[i][1])

print(c1)
print("c1,c2,b1,b2 are separated into distinct lists ")

parts = [[c1,c2],[b1,c1],[c2,b2],[c2,b1], [c1,b2], [b1,b2]]

strings = ["c1 and c2","b1 and c1","c2 and b2","c2 and b1","c1 and b2","b1 and b2"]

def get_combo(p1,p2):
    combo = []
    for i in range(len(c1)):
      combo.append(get_deltaR(p1[i],p2[i]))
      print("delta R "+str(get_deltaR(p1[i],p2[i])))
    combos.append(combo)

#computes minimum value
def get_min(hist_data):
    mean = np.mean(hist_data)
    std_dev = np.std(hist_data, ddof=1) 
    min = (mean - 2 * std_dev)-20
    return min 

#computes max value
def get_max(hist_data):
    mean = np.mean(hist_data)
    std_dev = np.std(hist_data, ddof=1) 
    min = (mean + 2 * std_dev)+20
    return min 

def make_hist(data, string, color):
    minimum = get_min(data)
    maximum = get_max(data)
    bins = 80
    
    hist = ROOT.TH1F(string+" hist", "delta between "+string, bins, minimum, maximum)
    for d in data:
        print(d)
        hist.Fill(d)

    hist.Scale(1.0 / hist.Integral())  
    hist.SetLineWidth(2)
    hist.SetLineColor(color)
    hist.SetMaximum(hist.GetMaximum() * 1.1) 
    hist.SetStats(0) 
    hist.Sumw2(ROOT.kFALSE)
    return hist

def make_plot(hists):
    canvas = ROOT.TCanvas("canvas", "deltaR between quark-anti-quark pairs", 800, 600)
    legend = ROOT.TLegend(0.1,0.7,0.48,0.9)

    for i, hist in enumerate(hists):
        mean = strings[i] + " deltaR mean: " + str(round(hist.GetMean(),2))
        legend.AddEntry(hist,mean,"l")

        if i ==0:
            hist.Draw("HIST")
        else:
            hist.Draw("SAME")
        
    legend.SetTextSize(0.0195)
    legend.SetBorderSize(0)
    legend.Draw()
    canvas.SaveAs("../finalHists/deltaR.pdf")

# contains lists of deltaRs for combinations
combos = []
for part in parts:
    get_combo(part[0],part[1])
    print("combos data is working")

hists = []
for i,combo in enumerate(combos):
    print(strings[i])
    print("does it get here?")
    hists.append(make_hist(combo, strings[i],colors[i]))
    print("or ?")

make_plot(hists)
