import uproot
import ROOT
import numpy as np
import awkward as ak 

rad = "1.0"
#root_file = "04corr30000.root"
#root_file = "chunk_1.root"
root_file = "1ecorr.root"


vars = ["jet_e_corr","recojet_e"]
plot_title = "Jet energy before vs. after c.o.m correction"
colors = [ROOT.kAzure+6,ROOT.kPink+9]

captions = ["after","before"]
types = ["total","higgs","zed","total summed"]

file = uproot.open(root_file)

print("finished reading")
    
tree = file['events']
branches = tree.arrays()

def get_var_array(var_name):
      return branches[var_name]
      

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

def get_total_event_e(jet_kind):
    jets = get_var_array(jet_kind)
    jets_mask = get_var_array("event_njet")==4
    mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
    mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)
    total_e = jets[jets_mask][mask]
    print("len of total energy array:",len(total_e))
    print_array(total_e)
    return jets[jets_mask][mask]
    
get_total_event_e(vars[1])