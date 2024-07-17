import ROOT 
import numpy as np
#read ROOT file
file = ROOT.TFile("ccH_Hbb_410_10k.root", "READ")
#get tree in Root File
antikt_tree = file.Get("events")

flavours = ["B","C","G","S","TAU","U"]
flavours_var = ["b","c","g","s","t","u"]

nbins_m = 50
xmin_m = 0.0
xmax_m = 25.0

nbins_p = 40
xmin_p = 0.0
xmax_p = 120.0

#create histogram from TH1 class
hist_m = ROOT.TH1F("hist_m", "Jet Mass; Jet Mass; Entries", nbins_m, xmin_m, xmax_m)
hist_p =  ROOT.TH1F("hist_p", "Jet Momentum; Jet Mass; Entries", nbins_p, xmin_p, xmax_p)
#iterate through the events in the three 
for event in antikt_tree:
    for i in range(len(antikt_tree.jet_mass)):
        m = event.jet_mass[i]
        p = event.jet_p[i]
        hist_p.Fill(p)
        hist_m.Fill(m)
 

#def pair_finder():


#events4 = []
    




#   #two_pairs=False
    #for a given event with 4 jets
    #flav1 is the flavor of first jet
    
    #flav1=ev.jets_truth[0]
    # #create vector for first jet flavour
    # f1=[]
    # #add first jet flavor to f1 vector
    # f1.append(flav1)
    # #iterate through the remaining jets
    # for i in range(1,len(ev)):
    #     #if jet flavor is equal to first jet flavor 
    #     if (ev.jets_truth[i]==flav1):
    #         #then add to f1 array 
    #         f1.append(ev.jets_truth[i])
    
    
    #if the f1 has a pair (exactly 2 elements)
    # if (nf1==2):
    #     if len(np.where(ev.jets_truth==):

    #     np.unique(0)
    #     #then for index i of elements in f1
    #     for f1[i] in f1:
    #         #if the index is not the 1
    #         if (i!=0):
    #             #save the index of the second element
    #             idx = i
    #         #initialize vector for remaining jets
    #     jets_remain =[]
    #         #for all the jets in the event
    #     for i in range(len(ev)):
    #         #if the jet is not the first flavour or second flavor
    #         if (i!=0 & i!=idx):
    #             #add to array of remaining jets
    #             jets_remain.append(ev.jets_truth[i])
    #     #if the flavor of remaining jets is equal
    #     if (jets_remain[0]==jets_remain[1]):
    #         #then there are 2 pairs with flav2 as flavour
    #         flav2 = jets_remain[0]
    #         two_pairs=True

    #     if (two_pairs) == True:
    #         if flav1 == 

        # for i in range(len(ev)):
        #     flavours=[]
        #     flavours.append(ev.jet_truth[i])

    #for i in range(len(events4)):
    
   
    # for i in range(len(antikt_tree.jets_truth)):
    #     t = event.jet_truth[i]
    #     for j in range(len(flavours)):
    #        flavours_var[i]= event.recojet_is[i]
        

##for each event-
    #if the number
##for each jet, jets_truth tells us what the pdgId of the flavour
##so for each event-- store all the jets_truth values in an array,
#hen 




#create canvas class 
canvas = ROOT.TCanvas("canvas", "Jet p", 800, 600)
hist_m.Draw()
canvas.SaveAs("histmp.pdf(")


canvas = ROOT.TCanvas("canvas", "Jet p", 800, 600)
hist_p.Draw()
canvas.SaveAs("histmp.pdf)")

