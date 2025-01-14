import uproot
import ROOT
import numpy as np
import awkward as ak
from matplotlib import pyplot as plt


files = ["1jconst.root", "1noreco.root","07ereco.root", "04reco.root"]
file_n=2
file = uproot.open(files[file_n])
tree = file['events']
branches = tree.arrays()

colors = ['red', 'magenta', 'blue', 'green']  
algs = ["anti-kt R=", "Durham-kt 4-jets mode"]
rads = ["0.4 ", "0.7 ", "1.0 "," "]
erecos = ["with energy recovery", "without energy recovery"," "]

alg = 0
rad = 1
ereco = 0

event_njet = branches["event_njet"]
print(event_njet)
print("length of event_njet prior to 4jet mask:", len(event_njet))

jets_mask = event_njet==4
print("jets mask:", jets_mask)
print("length jets mask:", len(jets_mask))
def get_event_indices(array):
    #initial index
    arr = np.array(array)
    event_idx=[]
    for i,ar in enumerate(arr):
        print(str(ar))
        if ar == True:
            event_idx.append(i)
    return event_idx

jet4events = branches["jets_truth"][jets_mask]
mask = np.array(ak.count_nonzero(np.abs(jet4events)==4, axis=1)==2)
mask &= np.array(ak.count_nonzero(np.abs(jet4events)==5, axis=1)==2)

events_index = get_event_indices(mask)
print("length of event index", len(events_index))
print("events index mask:", events_index)

print("mask length:", len(mask))
print("mask: ",mask)


b2c2events = jet4events[mask]
# print(len(get_event_indices(b2c2events)))
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
max_size=1000
min_size=0.0001

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

#define figure
fig = plt.figure(figsize=(10, 6), tight_layout=True, dpi=200)
ax = plt.subplot(1, 1, 1, aspect=1, xlim=[-np.pi, np.pi], ylim=[0, 3.8])

def reset_parameters(figure, axs):
    # Set custom ticks and labels
    axs.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    axs.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "$0$", r"$+\pi/2$", r"$+\pi$"])
    axs.set_yticks([0, np.pi / 2, np.pi])
    axs.set_yticklabels(["$0$", r"$+\pi/2$", r"$+\pi$"])
    axs.set_xlabel(r"$\mathbf{\phi}$")
    axs.set_ylabel(r"$\mathbf{\theta}$")

    # Customize spines
    axs.spines.right.set_visible(False)
    axs.spines.top.set_visible(False)
    axs.spines.left.set_position(("data", -np.pi - 0.2))
    axs.spines.left.set_bounds(0, np.pi)
    axs.spines.bottom.set_position(("data", -0.2))

reset_parameters(fig, ax)

#define scatter function
def get_scatter(phi, theta, markersize, color, transp, label):
    phi = np.array(phi)
    theta = np.array(theta)
    markersize = np.array(markersize)
    colors1= []
    for i in range(len(phi)):
        colors1.append(color)
    return ax.scatter(phi,theta,s=markersize, c=colors1, alpha=transp,label=label)

#define scatters function -- adds scatter plots to an array with their labels
def get_scatters(ind,labels):
    scatters = []
    for j in range(4):
        scatters.append(get_scatter(jc_phi_data[j][ind],jc_theta_data[j][ind], marker_sizes[j][ind], colors[j], transp,labels[j]))
        # print("scatters array made: event"+str(ind))
    return scatters

#define legend function
def get_legend(scatters): 
    return ax.legend(handles=scatters[:4], loc="upper right", title="Jet flavor", ncol=1, markerscale=0.4)


quarks = ["b","b","c","c"]

#point transparency
transp = 0.3


for i in range(dlength):
    get_legend(get_scatters(i,quarks))
    algo = algs[alg]
    radi = rads[rad]
    wreco = erecos[ereco]
    eventn = str(events_index[i]+1)
    event_number = str(i+1)
    ax.set_title("H(bb)Z(cc) event "+eventn+", "+algo+radi+wreco, loc='left',
                fontdict={'fontweight':'bold'})
                # , y=1.03)
    plt.show()
    plt.savefig("/usatlas/u/aconnelly/IzaFCCAnalysis/zhanalysis/plots/eventdisplays/antikt/ereco/07event"+event_number+
    ".png")
    print("plot"+event_number+
    ".png made")
    fig = plt.figure(figsize=(10, 6), tight_layout=True, dpi=200)
    ax = plt.subplot(1, 1, 1, aspect=1, xlim=[-np.pi, np.pi], ylim=[0, 3.8])
    reset_parameters(fig, ax)