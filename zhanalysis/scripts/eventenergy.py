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
    
def get_energies(energies_kind):
    event_energy = []
    print_array(energies_kind)
    for energies in energies_kind:
        event_energy.append(energies[0]+energies[1]+energies[2]+energies[3])
    print_array(event_energy)
    return event_energy
get_total_event_e(vars[0])

 
def make_hist(data,i,kind,max,min):
    maximum = max
    minimum = min
    label = types[kind]
    bins = 40
    hist = ROOT.TH1F("hist"+str(i),"zoom"+rad + " event " + label + " energy before vs. after c.o.m correction" , bins, minimum, maximum)
    for d in data:
        hist.Fill(d)
    hist.Scale(1.0 / hist.Integral())  
    hist.SetLineWidth(2)
    hist.SetLineColor(colors[i])
    hist.SetMaximum(hist.GetMaximum() * 1.1) 
    hist.SetStats(0) 
    hist.Sumw2(ROOT.kFALSE)
    return hist



def total_energy_hists(sums_list,kind,maxi,mini):
    label = types[kind]
    hists = []
    canvas = ROOT.TCanvas("canvas", "Event " + label + " energy before vs. after c.o.m correction", 1000, 1000)
    legend = ROOT.TLegend(0.1,0.7,0.48,0.9)
    

    for i,sum_list in enumerate(sums_list):
        hists.append(make_hist(sum_list,i,kind,maxi,mini))

    print(len(hists))

    for i, hist in enumerate(hists):
        mean = "energy " + captions[i] + " correction mean: " + str(round(hist.GetMean(),2))
        print(label+captions[i],hist.GetMean(),2)
        legend.AddEntry(hist,mean,"l")

        if i ==0:
            hist.Draw("HIST")
        else:
            hist.Draw("SAME")
        
    legend.SetTextSize(0.0195)
    legend.SetBorderSize(0)
    legend.Draw()

    canvas.SaveAs("../hists/"+"zoom"+str(rad) +"compare_"+label+"energies"+".pdf")


sums = [get_energies(get_total_event_e(vars[0])),get_energies(get_total_event_e(vars[1]))]

total_energy_hists(sums,0,270,160)

