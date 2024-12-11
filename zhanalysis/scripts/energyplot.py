import uproot
import ROOT
import numpy as np
import awkward as ak 

# variables = ["event_njet","jets_truth"]
# jet_corr= ["jet_e_corr","jet_px_corr","jet_py_corr","jet_pz_corr"]
# reco_jet = ["recojet_e","recojet_px","recojet_py","recojet_pz"]

rad = "dur"
# root_file = "04corr30000.root"
root_file = "chunk_1.root"

# before = "recojet_e"
# after = "jet_e_corr"

vars = ["jet_e_corr","recojet_e"]
plot_title = "Jet energy before vs. after c.o.m correction"
colors = [ROOT.kAzure+6,ROOT.kPink+9]

captions = ["after","before"]
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

#get array with all energies before and all energies after

def get_total_event_e(jet_kind):
    jets = get_var_array(jet_kind)
    jets_mask = get_var_array("event_njet")==4
    mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
    mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)
    return jets[jets_mask][mask]
    

def get_energies(energies_kind):
    event_energy = []
    for energies in energies_kind:
        event_energy.append(sum(energies))
    return event_energy


def get_energy_sums(energy_kind):
    e_when = get_total_event_e(energy_kind)
    e_sums_when =  get_energies(e_when)
    return e_sums_when


def make_hist(data,i):
    maximum = 300
    minimum = 0
    bins = 80
    hist = ROOT.TH1F("hist"+str(i), plot_title , bins, minimum, maximum)
    for d in data:
        hist.Fill(d)
    hist.Scale(1.0 / hist.Integral())  
    hist.SetLineWidth(2)
    hist.SetLineColor(colors[i])
    hist.SetMaximum(hist.GetMaximum() * 1.1) 
    hist.SetStats(0) 
    hist.Sumw2(ROOT.kFALSE)
    return hist
    

def total_energy_hists(sums_list):
    hists = []
    canvas = ROOT.TCanvas("canvas", "Jet energy before vs. after c.o.m correction", 1000, 1000)
    legend = ROOT.TLegend(0.1,0.7,0.48,0.9)

    for i,sum_list in enumerate(sums_list):
        hists.append(make_hist(sum_list,i))

    print(len(hists))

    for i, hist in enumerate(hists):
        mean = "energy " + captions[i] + " correction mean: " + str(round(hist.GetMean(),2))
        legend.AddEntry(hist,mean,"l")

        if i ==0:
            hist.Draw("HIST")
        else:
            hist.Draw("SAME")
        
    legend.SetTextSize(0.0195)
    legend.SetBorderSize(0)
    legend.Draw()

    canvas.SaveAs("../hists/"+str(rad) +"compare_energies"+".pdf")
    
sums_list=[]
for var in vars:
    sums_list.append(get_energy_sums(var))

total_energy_hists(sums_list)


