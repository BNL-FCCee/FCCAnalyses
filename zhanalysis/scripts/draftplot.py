import uproot
import ROOT
import numpy as np
from matplotlib import pyplot as plt


file0 = "1jconst.root"

caption = "Higgs truth energy from TLorentzVector"
label = "higgs_e_pepm_tlv"
particles = ["H","C","B"]
color = [ROOT.kRed-5]

file = uproot.open(file0)
tree = file['events']

branches = tree.arrays()
print("branches created")
colors= [ROOT.kRed-5,ROOT.kAzure+6, ROOT.kGreen+2]