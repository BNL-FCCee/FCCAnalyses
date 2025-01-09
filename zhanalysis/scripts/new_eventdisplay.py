import uproot
import ROOT
import numpy as np
import awkward as ak
from matplotlib import pyplot as plt

file0 = "1jconst.root"
file = uproot.open(file0)
tree = file['events']
branches = tree.arrays()
colors = ['red', 'magenta', 'blue', 'green']  

event_njet = branches["event_njet"]
print("length of event_njet prior to 4jet mask:", len(event_njet))

jets_mask = event_njet==4

jet4events = branches["jets_truth"][jets_mask]
print("length of event_njet after to 4jet mask:", len(jet4events))

mask = np.array(ak.count_nonzero(np.abs(jet4events)==4, axis=1)==2)
mask &= np.array(ak.count_nonzero(np.abs(jet4events)==5, axis=1)==2)

b2c2events = jet4events[mask]
print("length of event_njet after 2b and 2c:",len(b2c2events))
mask_c = np.abs(b2c2events)==4
mask_b = np.abs(b2c2events)==5

#jet constituents for jets from H
jc_phi_H = branches["jc_phi"][jets_mask][mask][mask_b]
jc_theta_H = branches["jc_theta"][jets_mask][mask][mask_b]
jc_e_H = branches["jc_e"][jets_mask][mask][mask_b]

#jet constituents for jets from Z
jc_phi_Z = branches["jc_phi"][jets_mask][mask][mask_c]
jc_theta_Z = branches["jc_theta"][jets_mask][mask][mask_c]
jc_e_Z = branches["jc_e"][jets_mask][mask][mask_c]

dlength = len(jc_phi_H)

jc_theta_data = []
jc_phi_data = []
jc_e_data = []

#each list in jc_data contains 4 arrays of jet constituents
jc_data = [jc_phi_data,jc_theta_data,jc_e_data]

#each H and Z lsit contain 2 arrays of jet constituents 
H_data = [jc_phi_H,jc_theta_H,jc_e_H]
Z_data = [jc_phi_Z,jc_theta_Z,jc_e_Z]

def get_particle_data(array, data_array):
    data1=[]
    data2=[]
    for i in range(len(array)):
        data1.append(array[i][0])
        data2.append(array[i][1])
    data_array.append(data1)
    data_array.append(data2)

#1st and 2nd data entries are b1 and b2 jet constituents
for i,data in enumerate(jc_data):
    get_particle_data(H_data[i],data)

#3rd and 4th data entries are c1 and c2 jet constituents
for i,data in enumerate(jc_data):
    get_particle_data(Z_data[i],data)

#each list in jc_data contains an list of lists containing b1 b2 c1 and c2 jet constituent data

marker_sizes = []
scale = 13
max_size=700
min_size=10

def get_marker_sizes(jc_energy):
    sizes = []
    for jc in jc_energy: 
        jc_size =[]
        for j in jc:
            size = np.pi * (np.sqrt(j)**2)*scale
            jc_size.append(size)
        jc_size1 = np.clip(jc_size, min_size, max_size)
        sizes.append(np.array(jc_size1))
    return sizes

#make marker sizes for b1, b2, c1, c2 data sets 
for jc_e in jc_e_data:
    marker_sizes.append(get_marker_sizes(jc_e))
    
leng = 20
def print_array(array):
    i=0
    # print(str(array))
    for arr in array:
        if i==leng:
            print(arr)
            # print(i)
            # print("length of arr array:", len(arr))
            print("done")
        if i<leng:
            print(arr)
            # print(i)
            # print("length of arr array:", len(arr))

        i = i+1

#define scatter function
def get_scatter(phi, theta, markersize, color, transp, axs):
    phi = np.array(phi)
    theta = np.array(theta)
    markersize = np.array(markersize)
    colors1= []
    for i in range(len(phi)):
        colors1.append(color)
    return axs.scatter(phi,theta,s=markersize, c=colors1, alpha=transp)

#define figure
fig = plt.figure(figsize=(10, 6), tight_layout=True, dpi=200)
ax = plt.subplot(1, 1, 1, aspect=1, xlim=[-np.pi, np.pi], ylim=[0, 3.8])

#point transparency
transp = 0.3

# Set custom ticks and labels
ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "$0$", r"$+\pi/2$", r"$+\pi$"])
ax.set_yticks([0, np.pi / 2, np.pi])
ax.set_yticklabels(["$0$", r"$+\pi/2$", r"$+\pi$"])
ax.set_xlabel(r"$\mathbf{\phi}$")
ax.set_ylabel(r"$\mathbf{\theta}$")

# Customize spines
ax.spines.right.set_visible(False)
ax.spines.top.set_visible(False)
ax.spines.left.set_position(("data", -np.pi - 0.2))
ax.spines.left.set_bounds(0, np.pi)
ax.spines.bottom.set_position(("data", -0.2))


# get_scatter(jc_phi_data[1][0],jc_theta_data[1][0], marker_sizes[1][0], colors[0], transp, ax)
# get_scatter(jc_phi_data[1][1],jc_theta_data[1][1], marker_sizes[1][1], colors[1], transp, ax)
# get_scatter(jc_phi_data[1][2],jc_theta_data[1][2], marker_sizes[1][2], colors[2], transp, ax)
# get_scatter(jc_phi_data[1][3],jc_theta_data[1][3], marker_sizes[1][3], colors[3], transp, ax)


get_scatter(jc_phi_data[0][2],jc_theta_data[0][2], marker_sizes[0][2], colors[0], transp, ax)

get_scatter(jc_phi_data[1][2],jc_theta_data[1][2], marker_sizes[1][2], colors[1], transp, ax)

get_scatter(jc_phi_data[2][2],jc_theta_data[2][2], marker_sizes[2][2], colors[2], transp, ax)

get_scatter(jc_phi_data[3][2],jc_theta_data[3][2], marker_sizes[3][2], colors[3], transp, ax)


# get_scatter(jc_phi_data[4][0],jc_theta_data[4][0], marker_sizes[4][0], colors[0], transp, ax)
# get_scatter(jc_phi_data[4][1],jc_theta_data[4][1], marker_sizes[4][1], colors[1], transp, ax)
# get_scatter(jc_phi_data[4][2],jc_theta_data[4][2], marker_sizes[4][2], colors[2], transp, ax)
# get_scatter(jc_phi_data[4][3],jc_theta_data[4][3], marker_sizes[4][3], colors[3], transp, ax)


plt.show()
plt.savefig("/usatlas/u/aconnelly/IzaFCCAnalysis/zhanalysis/plots/test8.png")
plt.close()

print("phi data lenngth",len(jc_phi_data[1]))
print("marker size data lenngth",len(marker_sizes[1]))



def print_array_length(array):
    i=0
    # print(str(array))
    for arr in array:
        if i==leng:
            for a in arr:
                print("length of a array:", len(a))
                print("a:",a)
            print("arr:", arr)
            print(i)
            print("done")
        if i<leng:
            for a in arr:
                print("length of a array:", len(a))
                print("a:",a)
            print("arr:", arr)
            print(i)

        i = i+1
