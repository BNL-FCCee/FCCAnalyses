import uproot
import ROOT
import numpy as np
import awkward as ak 

# event_start = 0
# #set number of events
# nevents = 4600
# #set number of bins
# bins = 40
# #legends for flavors
# legends =[[0.49, 0.7, 0.79, 0.9],[0.1,0.7,0.48,0.9]]
# compare = ["before", "after"]

variables = ["event_njet","jets_truth"]
jet_corr= ["jet_e_corr","jet_px_corr","jet_py_corr","jet_pz_corr"]
reco_jet = ["recojet_e","recojet_px","recojet_py","recojet_pz"]
root_file = "065corr.root"

before = "recojet_e"
after = "jet_e_corr"

captions = ["energy before","energy after"]
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

def get_event_energies(jet_kind):
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
    e_when = get_event_energies(energy_kind)
    e_sums_when =  get_energies(e_when)
    return e_sums_when


# print(get_energy_sums(before))
# print(get_energy_sums(after))

def plot_maker(sums_before, sums_after):
    # hists = []
    canvas = ROOT.TCanvas("canvas", "Jet energy before vs. after c.o.m correction", 1000, 1000)
    legend = ROOT.TLegend(0.49, 0.7, 0.79, 0.9)

    maximum = 300
    minimum = 150
    bins = 40

    # hists.append(hist)

    hist_a = ROOT.TH1F("hist_a","Jet energy before vs. after c.o.m correction", bins, minimum, maximum)
    for a in sums_after: 
        hist_a.Fill(a)
    hist_a.Scale(1.0 / hist_a.Integral()) 
    hist_a.SetLineColor(ROOT.kPink+9)
    hist_a.SetLineWidth(2)
    hist_a.SetMaximum(hist_a.GetMaximum() * 1.1) 
    hist_a.SetStats(0)
    hist_a.Sumw2(ROOT.kFALSE)
    mean_a = "Mean energy after:"+ str(round(hist_a.GetMean(),2))
    
    hist_a.Draw("HIST")


    hist = ROOT.TH1F("hist", "Jet energy before vs. after c.o.m correction" , bins, minimum, maximum)
    for b in sums_before: 
        hist.Fill(b)
    hist.Scale(1.0 / hist.Integral())  
    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kAzure+6)
    hist.SetMaximum(hist.GetMaximum() * 1.1) 
    hist.SetStats(0) 
    hist.Sumw2(ROOT.kFALSE)
    mean_b = "Mean before : "+str(round(hist.GetMean(),2))
    
    hist.Draw("SAME")

  
    # hists.append(hist_a)
    legend.AddEntry(hist,mean_b,"l")
    legend.AddEntry(hist_a,mean_a,"l")
    legend.SetTextSize(0.0195)
    legend.SetBorderSize(0)
    legend.Draw()

    canvas.SaveAs("../hists/"+"compare_energies"+".pdf")

    
plot_maker(get_energy_sums(before),get_energy_sums(after))

# e_before = get_event_energies("recojet_e")
# e_after = get_event_energies("jet_e_corr")
# e_sums_before = get_energies(e_before)
# e_sums_after = get_energies(e_after)


# print("energy before")
# print_array(e_sums_before)



