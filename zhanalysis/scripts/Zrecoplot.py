import uproot
import ROOT
import numpy as np
import awkward as ak 

#srad = same radius
#scorr = same corrections

# histogramStorage = ROOT.TFile.Open("histograms.root", "RECREATE")
# histogramStorage.Close()

#files = ["chunk_1.root","04corr30000.root", "06corr.root","06ptcorr.root","065corr.root", "07corr.root","09ecorr.root","1ecorr.root"]

#files = ["chunk_1.root","04corr30000.root","1ecorr.root"]
#files = ["chunk_1.root","04corr30000.root"]
# files = ["chunk_1.root","06corr.root","1ecorr.root"]
files = ["1ecorr.root"]



colors= [ROOT.kBlack, ROOT.kViolet+2,ROOT.kAzure+6]

flavs = ["cc","bb"]

# # key=str(rad[0])
       
type = [0,0,1]

event_start = 0
#set number of events
nevents = 10000
#set number of bins
bins = 40
#legends for flavors
legends =[[0.49, 0.7, 0.79, 0.9],[0.1,0.7,0.48,0.9]]

vars = ["px","py","pz","e"]
types = [["jet_","_corr"],"recojet_"]


# def srad_labels_maker(radius):
#     rad = str(radius)
#     labels = []
#     labels.append(" Durham 4jets with corr")
 
#     if type == 0:
#         labels.append(" R="+rad+" with ER and corr")
#     else:
#         labels.append(" R="+rad+" with ER")

#     return labels


def get_vars_p(branches,jets_mask,mask,mask_p):
    vs = []
    for var in vars:
        # if type == 0:
        #     v = branches[types[type][0]+var+types[type][1]]
        # if type == 1:
        #     v = branches[types[type]+var]
        v=branches["jet_"+var+"_corr"]
        va = v[jets_mask][mask][mask_p]
        vs.append(va)
    masses = get_masses(vs[0],vs[1],vs[2],vs[3])
    return masses


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

#idx is index of file
#first make all hists for each file

def make_hist(hhists,h_masses, idx):
            maximum = 200
            minimum = 20
            masses = h_masses
            print("length masses",len(masses))
            print(h_masses)
            hist = ROOT.TH1F("hist"+str(idx), "Higgs Masses", bins, minimum, maximum)
            for mass in masses: 
                for m in mass:
                    hist.Fill(m)
            hist.Scale(1.0 / hist.Integral())
            hhists.append(hist)



def create_plot(hists,flav,keys,lines):
    canvas = ROOT.TCanvas("canvas","Higgs Masses" , 800, 800)
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
        mean = str(round(hist.GetMean(),2))
        print(mean)
        legend.AddEntry(hist,keys[i]+"mean: "+mean,"l")
        legend.SetTextSize(0.026)
            
    legend.SetBorderSize(0)
    legend.Draw()
    canvas.SaveAs("../hists/rad_compare/FCCNote/Hmasses461.pdf")


def main(fi):    
    # lines = [2,1,1]
    lines = [1,1,1]

    # keys = ["Durham-kt 4-jet mode ", "Anti-kt with Jet Radius 0.6", "Anti-kt with Jet Radius 1.0"]
    keys = ["Anti-kt with Jet Radius 1.0 "]

    # cchists= []
    # bbhists= []
    # pphists=[cchists,bbhists]
    zmasses = []
    zhists = []
    for i,f in enumerate(fi):
        print(str(f))
        file = uproot.open(f)
        tree = file['events']
       
        print("finished reading")

        branches = tree.arrays()

        event_njet = branches["event_njet"]

        jets_mask = event_njet==4

        mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
        mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)
        mask_c = np.abs(branches["jets_truth"][jets_mask][mask])==4
        #mask_b = np.abs(branches["jets_truth"][jets_mask][mask])==5 

        #add Zed masses to hists
        zmasses.append(get_vars_p(branches,jets_mask,mask,mask_c))
        
        #added masses to collection of cc and bb masses
      
       # print(p_masses[0][:20])
        make_hist(zhists,zmasses,i)
        print("done with higgs")

    create_plot(zhists,1,keys,lines)
    



main(files)








