import uproot
import ROOT
import numpy as np
import awkward as ak

files = "r1truth.root"

file = uproot.open(files)
print("read file")
tree = file['events']

branches = tree.arrays(library = "np")
print("branches retrieved")

event_njet = branches["event_njet"]

jets_mask = event_njet==4
mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)
mask_c = np.abs(branches["jets_truth"][jets_mask][mask])==4
mask_b = np.abs(branches["jets_truth"][jets_mask][mask])==5 

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

def get_deltaR(data):
    for i in len(data[0]):
        p1 = ROOT.TLorentzVector()
        p2 = ROOT.TLorentzVector()
        p1.set(data[0][i][0],data[1][i][0],data[2][i][0],data[3][i][0])
        p2.set(data[0][i][1],data[1][i][1],data[2][i][1],data[3][i][1])
        

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

print_array(c_data[0])