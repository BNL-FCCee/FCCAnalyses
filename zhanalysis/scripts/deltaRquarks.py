import uproot
import ROOT
import numpy as np
import awkward as ak 

root_file = "1ecorrquark.root" 
file = uproot.open(root_file)
tree = file['events']
branches = tree.arrays()

##define function that takes eta and phi data from each quark
#then computes delta eta and delta phi and computes delta R
#save delta R in array 

def print_array(array):
    i=0
    for arr in array:
        if i==1000:
            print(arr)
            # print(i)
            print("done")
        if i<1000:
            print(arr)
            # print(i)
        i = i+1

truth_C_tlv = branches["truth_C_tlv"]
truth_C_e = branches["truth_C_e"]
truth_B_tlv = branches["truth_C_tlv"]
truth_B_e = branches["truth_B_e"]
truth_H_m = branches["truth_H_mass"]


print(truth_C_e)
print(truth_B_e)
print(truth_H_m)


# def get_phi(tlv):
#     phi = []
#     for tl in tlv:
#         phi_tl = []
#         for t in tl:
#             phi_tl.append(t.E())
#         phi.append(phi_tl)
#     return phi

phi_C = get_phi(truth_C_tlv)
