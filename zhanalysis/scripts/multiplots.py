import uproot
import ROOT
import numpy as np
import awkward as ak 

colors= [ROOT.kAzure+6, ROOT.kPink+9, ROOT.kGreen+2,
        ROOT.kOrange+3, ROOT.kMagenta+2,ROOT.kCyan+2,
        ROOT.kYellow+2, ROOT.kPink+1, ROOT.kBlack, 
        ROOT.kBlue-8]

flavs = ["cc","bb"]
#files = ["04nocorr.root","04corr30000.root", "06corr.root","065corr.root", "07corr.root","chunk_1.root"]
files = ["04corr30000.root", "06corr.root","065corr.root", "07corr.root","chunk_1.root"]

key = "allplots"

#set starting event
event_start = 0
#set number of events
nevents = 500
#set number of bins
bins = 35

legends =[[0.49, 0.7, 0.79, 0.9],[0.1,0.7,0.48,0.9]]

labels=[[" 0.4 with ER", " 0.4 with ER and corr"],
        [" 0.6 with ER"," 0.6 with ER and corr"],
        ["0.65 with ER"," 0.65 with ER and corr"],
        [ " 0.7 with ER", " 0.7 with ER and corr"],
        [" Durham 4jets with corr"]]


howmass = [2,1,2] 
lines = [1,1,2]

howmass = [2,2,2,2,2]
lines = [1,1,1,1,2]

howmass = [1,1,1,1,2]

labels_c ={}
labels_b = {}
labels_p = [labels_c,labels_b]

cchists= []
bbhists= []
pphists=[cchists,bbhists]


def get_masses(px,py,pz,e):
    p_mass=[]
    for i in range(len(px)):
    #for i in range(nevents):
        vecs=[]

        vecs.append(ROOT.TLorentzVector())
        vecs.append(ROOT.TLorentzVector())

        vecs[0].SetPxPyPzE(px[i][0],py[i][0], pz[i][0], e[i][0])
        vecs[1].SetPxPyPzE(px[i][1],py[i][1], pz[i][1], e[i][1])

        sum = ROOT.TLorentzVector()
        #compute tlv sum
        sum = vecs[0] + vecs[1]
    
        #get invariant mass of the sum 
        mass=sum.M()

        p_mass.append(mass)
        if i % 10000 == 0:
         print(i)
        
    return p_mass

def assign_vars(p_masses, howmass, branches, jets_mask, mask, mask_c, mask_b):
   
   if howmass==0:
    jet_px_c=(branches["recojet_px"][jets_mask][mask][mask_c])
    jet_py_c=(branches["recojet_py"][jets_mask][mask][mask_c])
    jet_pz_c=(branches["recojet_pz"][jets_mask][mask][mask_c])
    jet_e_c=(branches["recojet_e"][jets_mask][mask][mask_c])

    jet_px_b=(branches["recojet_px"][jets_mask][mask][mask_b])
    jet_py_b=(branches["recojet_py"][jets_mask][mask][mask_b])
    jet_pz_b=(branches["recojet_pz"][jets_mask][mask][mask_b])
    jet_e_b=(branches["recojet_e"][jets_mask][mask][mask_b])

    c_masses = get_masses(jet_px_c,jet_py_c,jet_pz_c,jet_e_c)
    b_masses = get_masses(jet_px_b,jet_py_b,jet_pz_b,jet_e_b)
   
    p_masses.append(c_masses)
    p_masses.append(b_masses)
   
   elif howmass==1:
    jet_px_ccorr=(branches["jet_px_corr"][jets_mask][mask][mask_c])
    jet_py_ccorr=(branches["jet_py_corr"][jets_mask][mask][mask_c])
    jet_pz_ccorr=(branches["jet_pz_corr"][jets_mask][mask][mask_c])
    jet_e_ccorr=(branches["jet_e_corr"][jets_mask][mask][mask_c])

    jet_px_c=(branches["recojet_px"][jets_mask][mask][mask_c])
    jet_py_c=(branches["recojet_py"][jets_mask][mask][mask_c])
    jet_pz_c=(branches["recojet_pz"][jets_mask][mask][mask_c])
    jet_e_c=(branches["recojet_e"][jets_mask][mask][mask_c])


    jet_px_bcorr=(branches["jet_px_corr"][jets_mask][mask][mask_b])
    jet_py_bcorr=(branches["jet_py_corr"][jets_mask][mask][mask_b])
    jet_pz_bcorr=(branches["jet_pz_corr"][jets_mask][mask][mask_b])
    jet_e_bcorr=(branches["jet_e_corr"][jets_mask][mask][mask_b])
   
    jet_px_b=(branches["recojet_px"][jets_mask][mask][mask_b])
    jet_py_b=(branches["recojet_py"][jets_mask][mask][mask_b])
    jet_pz_b=(branches["recojet_pz"][jets_mask][mask][mask_b])
    jet_e_b=(branches["recojet_e"][jets_mask][mask][mask_b])


    c_masses_corr = get_masses(jet_px_ccorr,jet_py_ccorr,jet_pz_ccorr,jet_e_ccorr)
    c_masses = get_masses(jet_px_c,jet_py_c,jet_pz_c,jet_e_c)

    b_masses_corr = get_masses(jet_px_bcorr,jet_py_bcorr,jet_pz_bcorr,jet_e_bcorr)
    b_masses = get_masses(jet_px_b,jet_py_b,jet_pz_b,jet_e_b)
    
    p_masses.append(c_masses)
    p_masses.append(c_masses_corr)

    p_masses.append(b_masses)
    p_masses.append(b_masses_corr)
    
    
   else: 
    jet_px_c=(branches["jet_px_corr"][jets_mask][mask][mask_c])
    jet_py_c=(branches["jet_py_corr"][jets_mask][mask][mask_c])
    jet_pz_c=(branches["jet_pz_corr"][jets_mask][mask][mask_c])
    jet_e_c=(branches["jet_e_corr"][jets_mask][mask][mask_c])

    jet_px_b=(branches["jet_px_corr"][jets_mask][mask][mask_b])
    jet_py_b=(branches["jet_py_corr"][jets_mask][mask][mask_b])
    jet_pz_b=(branches["jet_pz_corr"][jets_mask][mask][mask_b])
    jet_e_b=(branches["jet_e_corr"][jets_mask][mask][mask_b])

    c_masses = get_masses(jet_px_c,jet_py_c,jet_pz_c,jet_e_c)
    b_masses = get_masses(jet_px_b,jet_py_b,jet_pz_b,jet_e_b)

    p_masses.append(c_masses)
    p_masses.append(b_masses)

#idx is index of file
#first make all hists for the 
def add_file_hists(p_masses, flav, idx):
    for i in range(len(labels[idx])):
        maximum = 200
        minimum = 0
        if len(p_masses)==4 and flav==0:
            masses = p_masses[i]
        elif len(p_masses)==4 and flav==1:
            masses = p_masses[i+2]
        else:
            masses = p_masses[flav]

        hist = ROOT.TH1F("hist"+str(idx)+str(i)+flavs[flav],flavs[flav]+ " jet pair masses with anti-kt jet clustering", bins, minimum, maximum)
        for m in masses: 
            hist.Fill(m)
        mean = "Mean: "+str(round(hist.GetMean(),2))
        #means_p[flav][hist]=mean
        label = labels[idx][i]
        labels_p[flav][hist] =mean+label
        #color = colors[idx+i]
        hist.Scale(1.0 / hist.Integral())
        pphists[flav].append(hist)


def create_plot(hists,flav):
    canvas = ROOT.TCanvas("canvas", flavs[flav]+" jet pair masses", 1000, 1000)

    legend = ROOT.TLegend(legends[flav][0],legends[flav][1],legends[flav][2],legends[flav][3])
       
    
    for i, hist in enumerate(hists):
        color = colors[i]
        hist.SetLineColor(color)
        hist.SetLineWidth(2)
        hist.SetLineStyle(lines[i])
        hist.SetMaximum(hist.GetMaximum() * 1.1) 
        hist.SetStats(0) 
        hist.Sumw2(ROOT.kFALSE)
        if i==0:
            hist.Draw("HIST")
        else:
            hist.Draw("SAME")
        #legend.AddEntry(hist,labels[i]+ mean,"l")
        legend.SetTextSize(0.0195)
        legend.AddEntry(hist,labels_p[flav][hist],"l")
        
        #legend.AddEntry(hist,means_p[flav][hist],"p")
    
    legend.SetBorderSize(0)
    legend.Draw()

    canvas.SaveAs("../hists/"+key+flavs[flav]+"mass.pdf")


def main(f,idx):
    file = uproot.open(f)
    tree = file['events']
    #array of the branches
    p_masses=[] 
    print("i read it hehe")
    branches = tree.arrays()

    event_njet = branches["event_njet"]

    if any(event_njet)>4:
        print("something is wrong")

    jets_mask = event_njet==4

    mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
    mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)

    mask_c = np.abs(branches["jets_truth"][jets_mask][mask])==4
    mask_b = np.abs(branches["jets_truth"][jets_mask][mask])==5

    assign_vars(p_masses, howmass[idx], branches, jets_mask, mask, mask_c, mask_b)

    add_file_hists(p_masses,0,idx)
    add_file_hists(p_masses,1,idx)

  


for i,f in enumerate(files): 
    main(f,i)

for i,hists in enumerate(pphists): 
    create_plot(hists,i)

    




