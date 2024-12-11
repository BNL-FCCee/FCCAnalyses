import uproot
import ROOT
import numpy as np


files = "r1truth.root"

file = uproot.open(files)
print("read file")
tree = file['events']

branches = tree.arrays(library = "np")
print("branches retrieved")

colors= [ROOT.kRed-5,ROOT.kAzure+6, ROOT.kGreen+2]
filename= "higgs_pt"

def print_array(array):
    i=0
    for arr in array:
        if i==100:
            print(arr)
            # print(i)
            print("done")
        if i<100:
            print(arr)
            # print(i)
        i = i+1


# vars = ["p","pt","px","py","pz","phi","mass","eta"]
#computes minimum value
def get_min(hist_data):
    mean = np.mean(hist_data)
    std_dev = np.std(hist_data, ddof=1) 
    min = (mean - 2 * std_dev)-5
    return min 

#computes max value
def get_max(hist_data):
    mean = np.mean(hist_data)
    std_dev = np.std(hist_data, ddof=1) 
    min = (mean + 2 * std_dev)+5
    return min 

p = branches["truth_H_pt"]
# print(p[:10])
# m2p = []
# for i in range(100000):
#     data = (2*125)/(p[i][0])
#     m2p.append(data)


# for i in range(10):
#     print(p[i][0])

def make_hist(data):
    
    canvas = ROOT.TCanvas("canvas", "higgs vars", 800, 600)

    # minimum = get_min(data)
    minimum=0
    maximum = get_max(data)
    bins = 80
    
    hist = ROOT.TH1F("higgs var", "higgs pt", bins, minimum, maximum)

    for dat in data:
        # hist.Fill(dat)
        for d in dat:
            hist.Fill(d)

    hist.Scale(1.0 / hist.Integral())
    hist.SetMaximum(hist.GetMaximum() * 1.1) 
    hist.Sumw2(ROOT.kFALSE)
    color = colors[0]
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.Draw("HIST")
    canvas.SaveAs("../hists/higgsvars/"+filename+".pdf")

# make_hist(p)
make_hist(p)