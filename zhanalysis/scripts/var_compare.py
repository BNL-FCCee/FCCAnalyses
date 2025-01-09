import uproot
import ROOT
import numpy as np
import awkward as ak 

#srad = same radius
#scorr = same corrections

histogramStorage = ROOT.TFile.Open("histograms.root", "RECREATE")
histogramStorage.Close()

#files = ["chunk_1.root","04corr30000.root", "06corr.root","06ptcorr.root","065corr.root", "07corr.root","09ptcorr.root","1ecorr.root"]

#files = ["chunk_1.root","1ecorr.root"]

file = "hdecay.root"


colors= [ROOT.kBlack, ROOT.kViolet+2]

flavs = ["cc","bb"]

rad=0.4

key=" "+str(rad)+" "   

event_start = 0
#set number of events
nevents = 4600
#set number of bins
bins = 40
#legends for flavors
legends =[[0.49, 0.7, 0.79, 0.9],[0.1,0.7,0.48,0.9]]


vars = ["px","py","pz","e"]
names = ["recojet_",["jet_","_corr"]]
types = [0,1]
compare = ["before", "after"]
variables = ["mass ","energy "]

# def srad_labels_maker(radius):
#     rad = str(radius)
#     labels = []
#     labels.append(" Durham 4jets with corr")
 
#     if type == 0:
#         labels.append(" R="+rad+" with ER")
#     else:
#          labels.append(" R="+rad+" with ER and corr")

#     return labels

def get_masses(px,py,pz,e):
    p_mass=[]
    #for i in range(len(px)):
    for i in range(nevents):
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

def get_energies(e):
    energy_per_event = []
    print (len(e))
    for i in range(nevents):
        energy_sum = e[i][0] + e[i][1]+ e[i][2]+e[i][3]
        energy_per_event.append(energy_sum)
    print("energy per event length", len(energy_per_event)  )
    return energy_per_event
        

def get_e(branches,jets_mask,mask,kind):
    if kind == 0:
        e = branches[names[kind]+"e"]
        
    if kind == 1:
        e = branches[names[kind][0]+"e"+names[kind][1]]
   
    en = e[jets_mask][mask]

    e_energies = get_energies(en)
    return e_energies


def get_vars_masses(branches,jets_mask,mask,mask_p,kind):
    vs = []
    for var in vars:
        if kind == 0:
            v = branches[names[kind]+var]
            
        if kind == 1:
            v = branches[names[kind][0]+var+names[kind][1]]
        va = v[jets_mask][mask][mask_p]
        vs.append(va)
    v_masses = get_masses(*vs)
    return v_masses


#gets reco masses or corr masses for both flavors
def get_type_masses(branches,jets_mask,mask,masks_p,kind):
        type_masses =[]
        z_masses = get_vars_masses(branches, jets_mask,mask,masks_p[0],kind)
        h_masses = get_vars_masses(branches, jets_mask,mask,masks_p[1],kind)
        type_masses.append(z_masses)
        type_masses.append(h_masses)
        return type_masses

#idx is index of file
#first make all hists for each file

#def add_file_hists(pphists,p_masses, flav, key, idx):
#             maximum = 200
#             minimum = 0
#             masses = p_masses[flav]
#             print("length masses",len(masses))
#             for i,p_mass in enumerate(masses):
#                 hist = ROOT.TH1F("hist"+str(idx)+str(i)+flavs[flav], "R="+key , bins, minimum, maximum)
#                 for pp in p_mass: 
#                     hist.Fill(pp)
#                 #color = colors[idx+i]
#                 histogramStorage = ROOT.TFile.Open("histograms.root","UPDATE")
#                 histogramStorage.cd()
#                 hist.Write()
#                 histogramStorage.Close()
#                 hist.Scale(1.0 / hist.Integral())
#                 pphists[flav].append(hist)


#def create_plot(hists,flav,labels,key,lines):
#     canvas = ROOT.TCanvas("canvas", "radius "+key , 800, 800)
#     legend = ROOT.TLegend(legends[flav][0],legends[flav][1],legends[flav][2],legends[flav][3])
    
#     for i, hist in enumerate(hists):
#         color = colors[i]
#         hist.SetLineColor(color)
#         hist.SetLineWidth(2)
#         hist.SetLineStyle(lines[i])
#         hist.SetMaximum(hist.GetMaximum() * 1.1) 
#         hist.SetStats(0) 
#         hist.Sumw2(ROOT.kFALSE)
#         if i==0:
#             hist.Draw("HIST")
#         else:
#             hist.Draw("SAME")
#         mean = "Mean: "+str(round(hist.GetMean(),2))
#         print(mean)
#         legend.AddEntry(hist,mean+labels[i],"l")
#         legend.SetTextSize(0.026)
            
#     legend.SetBorderSize(0)
#     legend.Draw()
#     canvas.SaveAs("../hists/"+key+flavs[flav]+"mass.pdf")

def get_plot_vars(fi):   
   # cc masses before and after corr 
    reco_mass_hz=[]
    corr_mass_hz = [] 
  
    reco_energy = []
    corr_energy = []

    file = uproot.open(fi)
    print("finished reading")

    tree = file['events']
  
    branches = tree.arrays()

    event_njet = branches["event_njet"]

    jets_mask = event_njet==4
    mask = np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==4, axis=1)==2)
    mask &= np.array(ak.count_nonzero(np.abs(branches["jets_truth"][jets_mask])==5, axis=1)==2)

    mask_c = np.abs(branches["jets_truth"][jets_mask][mask])==4
    mask_b = np.abs(branches["jets_truth"][jets_mask][mask])==5

    masks_p = [mask_c,mask_b]

    reco_mass_hz.append(get_type_masses(branches,jets_mask,mask,masks_p,0))
    corr_mass_hz.append(get_type_masses(branches,jets_mask,mask,masks_p,1))

    reco_energy.append(get_e(branches,jets_mask,mask,0))
    corr_energy.append(get_e(branches,jets_mask,mask,1))

    # for type in types:
    #     type_hz[type].append(get_type_masses(branches,jets_mask,mask,masks_p,type))
    
    type_hz = [reco_mass_hz,corr_mass_hz,reco_energy, corr_energy]
    
    return type_hz

def make_scatter_plot(data,vari,when):
    c1 = ROOT.TCanvas( 'c1', 'comparison', 200, 10, 700, 500)
    # c1.SetGrid()
    x = np.array(data[0],dtype=np.float64)
    y = np.array(data[1],dtype=np.float64)
    #e.g. mass before energy correction, energy after correction
    label = variables[vari]+compare[when]+" correction"
    graph = ROOT.TGraph(nevents,x,y)
    # graph.SetLineColor(2)
    # graph.SetLineWidth(4)
    graph.SetMarkerColor(ROOT.kBlue-2)
    graph.SetMarkerStyle(7)
    graph.SetTitle(key+variables[vari]+compare[when])
    graph.GetXaxis().SetTitle("Z "+label)
    graph.GetYaxis().SetTitle("H "+label)
    graph.Draw("AP")
    c1.SaveAs("../plots/"+key+variables[vari]+compare[when]+".pdf")

hztype = get_plot_vars(file)

data_m = [hztype[0],hztype[1]]
data_e = [hztype[2],hztype[3]]

print(len(data_e[0]))
print(len(data_e[0][0]))

print(len(data_m[0]))
print(len(data_m[0][0]))
print(len(data_m[0][0][0]))


#0 for before and 1 for after
for i in range(2):
    make_scatter_plot(data_e,1,i)


for i in range(2):
    make_scatter_plot(data_m[i][0],0,i)






