import uproot
import ROOT
import numpy as np

files = ["1ecorrquark.root"]

caption = "Higgs truth energy from TLorentzVector"
label = "higgs_e_pepm_tlv"
particles = ["H","C","B"]
color = [ROOT.kRed-5]

file = uproot.open(files[0])
tree = file['events']

branches = tree.arrays()
colors= [ROOT.kRed-5,ROOT.kAzure+6, ROOT.kGreen+2]


def print_array(array):
    i=0
    for arr in array:
        if i==100:
            print(arr)
            print(i)
            print("done")
        if i<100:
            print(arr)
            print(i)
        i = i+1

tvars = ["pt","eta","phi","mass","e","phi","eta"]
t_vars=[]

def assign_branches(particle):
    t_vars=[]
    for var in tvars:
        t_vars.append(branches["truth_"+particle+"_"+var])
    return t_vars

for particle in particles:
    assign_branches(particle)

print("C phi")
print(branches["truth_C_phi"])
print("C eta")
print(branches["truth_C_eta"])
print("B eta")
print(branches["truth_B_eta"])
print("B phi")
print(branches["truth_B_phi"])
print("C energy")
print(branches["truth_C_e"])
print("B energy")
print(branches["truth_B_e"])
print("H phi")
print(branches["truth_H_phi"])

# print("t_vars 0 0 0:", t_vars[0][0][0])
# print("t_vars 0 1 0:", t_vars[0][1][0])

# print("length of t_vars",len(t_vars))
# print("length of t_vars[0]", len(t_vars[0]))

#get Higgs from TLV
def get_energies(pt, eta, phi, mass):
    e_tlv=[]
    for i in range(len(pt)):
        tlv = ROOT.TLorentzVector()
        tlv.SetPtEtaPhiM(pt[i][0], eta[i][0], phi[i][0], mass[i][0])
        e_tlv.append(tlv.E())
    return e_tlv

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


def make_hist(data):
    
    canvas = ROOT.TCanvas("canvas", "higgs energy from TLV", 800, 600)

    hist_data = data

    minimum = get_min(hist_data)
    maximum = get_max(hist_data)
    bins = 80
    
    hist = ROOT.TH1F("hist_var", "higgs energy from TLV", bins, minimum, maximum)

    for h_d in hist_data:
        hist.Fill(h_d)

    hist.Scale(1.0 / hist.Integral())
    hist.SetMaximum(hist.GetMaximum() * 1.1) 
    hist.Sumw2(ROOT.kFALSE)
    color = colors[0]
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.Draw("HIST")
    canvas.SaveAs("../hists/higgsvars/"+label+".pdf")

# higgs_e = get_energies(t_vars[0],t_vars[1],t_vars[2],t_vars[3])


# make_hist(branches["truth_"])

# make_hist(higgs_e)