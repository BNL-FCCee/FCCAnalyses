import uproot
import ROOT
import numpy as np

colors= [ROOT.kRed-5,ROOT.kAzure+6, ROOT.kGreen+2]

files = ["1ecorr.root"]

particle = "H"
file = uproot.open(files[0])
tree = file['events']

branches = tree.arrays()

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


vars = ["e","p","pt","px","py","pz","phi","mass","eta"]
tvars = ["pt","eta","phi","mass"]

p_vars=[]
for var in vars: 
    p_vars.append(branches["truth_"+particle+"_"+var])

t_vars=[]
for var in tvars:
    t_vars.append(branches["truth_"+particle+"_"+var])

print("higgs vars array format: ",p_vars)

#get Higgs from TLV
def get_energies(pt, eta, phi, mass):
    e_tlv=[]
    for i in len(t_vars[0]):
        tlv = ROOT.TLorentzVector
        tlv.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i])
        e_tlv.append(tlv.E())
    return e_tlv
        


def get_data(p_vars):
    hist_data =[]
    # print_array(higgs_vars)
    for val in p_vars:
        # print_array(val)
        for v in val:
            hist_data.append(v)
    return hist_data


#0 for min and 1 for max
def get_min(hist_data):
    mean = np.mean(hist_data)
    std_dev = np.std(hist_data, ddof=1) 
    min = (mean - 2 * std_dev)-10
    return min 

def get_max(hist_data):
    mean = np.mean(hist_data)
    std_dev = np.std(hist_data, ddof=1) 
    min = (mean + 2 * std_dev)+10
    return min 


def make_hist(p_vars, var):
    caption = "particle"+var

    canvas = ROOT.TCanvas("canvas", "higgs "+var, 800, 600)

    hist_data = get_data(p_vars)

    minimum = get_min(hist_data)
    maximum = get_max(hist_data)
    bins = 80
    
    hist = ROOT.TH1F("hist_var", "particle "+var, bins, minimum, maximum)

    for data in hist_data:
        hist.Fill(data)

    # for var in higgs_vars:
    #     for v in var:
    #         hist.Fill(v)

    hist.Scale(1.0 / hist.Integral())
    hist.SetMaximum(hist.GetMaximum() * 1.1) 
    hist.Sumw2(ROOT.kFALSE)
    color = colors[0]
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.Draw("HIST")
    canvas.SaveAs("../hists/higgsvars/"+"70bins_"+caption+".pdf")


# for i,p_var in enumerate(p_vars):
#     make_hist(p_var,vars[i])




#make_hist(higgs_vars[3],vars[3])

