import uproot
import ROOT
import numpy as np

colors= [ROOT.kAzure+6, ROOT.kPink+9, ROOT.kGreen+2,
        ROOT.kOrange+3, ROOT.kMagenta+2,ROOT.kCyan+2,
        ROOT.kYellow+2, ROOT.kPink+1, ROOT.kBlack, 
        ROOT.kBlue-8]

files = ["06corr.root"]

key = "06pt"

file = uproot.open(files[0])
tree = file['events']

branches = tree.arrays()

vars = ["p","m","e","pt","px","py","pz","phi","eta"]
higgs_vars=[]
higgs_vars.append(branches["truth_H_p"])
higgs_vars.append(branches["truth_H_mass"])
higgs_vars.append(branches["truth_H_e"])
higgs_vars.append(branches["truth_H_pt"])
higgs_vars.append(branches["truth_H_px"])
higgs_vars.append(branches["truth_H_py"])
higgs_vars.append(branches["truth_H_pz"])
higgs_vars.append(branches["truth_H_phi"])
higgs_vars.append(branches["truth_H_eta"])

overh = []

def make_hist(higgs_vars, vars):
    key = vars

    h_values =[]
    for val in higgs_vars:
        for v in val:
            h_values.append(v)
            
    maximum = 100
    minimum = 0
    bins = 40

    
    hist = ROOT.TH1F("hist_var", "higgs"+vars, bins, minimum, maximum)

    for var in higgs_vars:
        for v in var:
            hist.Fill(v)

    hist.Scale(1.0 / hist.Integral())
    hist.SetMaximum(hist.GetMaximum() * 1.1) 
    hist.Sumw2(ROOT.kFALSE)
    color = colors[1]
    hist.SetLineColor(color)
    hist.SetLineWidth(2)

    canvas = ROOT.TCanvas("canvas", "higgs"+vars, 800, 600)


    hist.Draw("HIST")

    canvas.SaveAs("../hists/"+key+".pdf")

#make_hist(higgs_vars[3],vars[3])

print(higgs_vars[0])