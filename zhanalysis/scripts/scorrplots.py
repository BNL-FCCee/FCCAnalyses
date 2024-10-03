import uproot
import ROOT
import numpy as np
import awkward as ak 

#srad = same radius
#scorr = same corrections

# files = ["chunk_1.root","04corr30000.root", "06corr.root","065corr.root", "07corr.root"]
files = ["chunk_1.root","04corr30000.root","1ecorr.root"]

colors= [ROOT.kAzure+6, ROOT.kPink+9, ROOT.kGreen+2,
        ROOT.kOrange+3, ROOT.kMagenta+2,ROOT.kCyan+2,
        ROOT.kYellow+2, ROOT.kPink+1, ROOT.kBlack, 
        ROOT.kBlue-8]

flavs = ["cc","bb"]

corrs = ["withERanddcorr","withER"]


#set starting event
event_start = 0
#set number of events
nevents = 500
#set number of bins
bins = 35

legends =[[0.49, 0.7, 0.79, 0.9],[0.1,0.7,0.48,0.9]]


kinds = ["with ER and corr","with ER"]

#0 -- energy recovery 1 -- energy recovery and correction
def scorr_labels_maker(kind):
    corre = kinds[kind]
    labels = ["Durham 4jets with dcorr"," 0.4"+corre," 0.6"+corre," 0.65"+corre, " 0.75"+corre]
    return labels


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


vars = ["px","py","pz","e"]
types = [["jet_","_corr"],"recojet_"]

def get_vars_p(branches,jets_mask,mask,mask_p,type):
    vs = []
    for var in vars:
        if type == 0:
            v = branches[types[type][0]+var+types[type][1]]
        if type == 1:
            v = branches[types[type]+var]
        va = v[jets_mask][mask][mask_p]
        vs.append(va)
    masses = get_masses(vs[0],vs[1],vs[2],vs[3])
    return masses


#idx is index of file
#first make all hists for each file


def make_plot(data):
    plot = ROOT.TGraph()
    


def add_file_hists(p_masses, flav, idx):
            maximum = 200
            minimum = 0
            masses = p_masses[flav]
            for i,p_mass in enumerate(masses):
                hist = ROOT.TH1F("hist"+str(idx)+str(i)+flavs[flav],flavs[flav]+ " jet pair masses with anti-kt jet clustering", bins, minimum, maximum)
                for p in p_mass: 
                    hist.Fill(p)
                #color = colors[idx+i]
                hist.Scale(1.0 / hist.Integral())
                p_masses[flav].append(hist)


def create_plot(hists,flav,labels,key,lines):
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
        mean = "Mean: "+str(round(hist.GetMean(),2))
        legend.AddEntry(hist,mean+labels[i],"l")
        legend.SetTextSize(0.0195)
            
    legend.SetBorderSize(0)
    legend.Draw()
    canvas.SaveAs("../hists/"+key+flavs[flav]+"mass.pdf")


def main(fi,corr):    
    lines = [2,1,1,1,1]
    cchists= []
    bbhists= []
    pphists=[cchists,bbhists]

    for i,f in enumerate(fi):
        p_masses=[]
        c_masses=[]
        b_masses=[]

        file = uproot.open(f)
        tree = file['events']

        print("finished reading")

        branches = tree.arrays()

        event_njet = branches["event_njet"]

        jets_mask = event_njet==4

        mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
        mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)
        mask_c = np.abs(branches["jets_truth"][jets_mask][mask])==4
        mask_b = np.abs(branches["jets_truth"][jets_mask][mask])==5

        #added masses to collection of cc and bb masses

        if corr == 0:
            key = corrs[corr]
            if i>=0:
                c_masses.append(get_vars_p(branches,jets_mask,mask,mask_c,0))
                b_masses.append(get_vars_p(branches,jets_mask,mask,mask_b,0))

        if corr == 1:
            key = corrs[corr]
            if i == 0:
                c_masses.append(get_vars_p(branches,jets_mask,mask,mask_c,0))
                b_masses.append(get_vars_p(branches,jets_mask,mask,mask_b,0))
            
            if i>0:
                c_masses.append(get_vars_p(branches,jets_mask,mask,mask_c,1))
                b_masses.append(get_vars_p(branches,jets_mask,mask,mask_b,1))

        p_masses.append(c_masses)
        p_masses.append(b_masses)

        print(p_masses[0][:20])

        add_file_hists(p_masses,0,i)
        add_file_hists(p_masses,1,i)

        labels = scorr_labels_maker(corr)

        for i,pphist in enumerate(pphists):
            create_plot(pphist,i,labels,key,lines)


main(files,0)
main(files,1)







