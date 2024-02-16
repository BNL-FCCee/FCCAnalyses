#!/usr/bin/env python

import ROOT
import numpy as np
from matplotlib import pyplot as plt

colors = {
    "gray": {
        0: "#f8f9fa",
        1: "#f1f3f5",
        2: "#e9ecef",
        3: "#dee2e6",
        4: "#ced4da",
        5: "#adb5bd",
        6: "#868e96",
        7: "#495057",
        8: "#343a40",
        9: "#212529",
    },
    "red": {
        0: "#fff5f5",
        1: "#ffe3e3",
        2: "#ffc9c9",
        3: "#ffa8a8",
        4: "#ff8787",
        5: "#ff6b6b",
        6: "#fa5252",
        7: "#f03e3e",
        8: "#e03131",
        9: "#c92a2a",
    },
    "pink": {
        0: "#fff0f6",
        1: "#ffdeeb",
        2: "#fcc2d7",
        3: "#faa2c1",
        4: "#f783ac",
        5: "#f06595",
        6: "#e64980",
        7: "#d6336c",
        8: "#c2255c",
        9: "#a61e4d",
    },
    "grape": {
        0: "#f8f0fc",
        1: "#f3d9fa",
        2: "#eebefa",
        3: "#e599f7",
        4: "#da77f2",
        5: "#cc5de8",
        6: "#be4bdb",
        7: "#ae3ec9",
        8: "#9c36b5",
        9: "#862e9c",
    },
    "violet": {
        0: "#f3f0ff",
        1: "#e5dbff",
        2: "#d0bfff",
        3: "#b197fc",
        4: "#9775fa",
        5: "#845ef7",
        6: "#7950f2",
        7: "#7048e8",
        8: "#6741d9",
        9: "#5f3dc4",
    },
    "indigo": {
        0: "#edf2ff",
        1: "#dbe4ff",
        2: "#bac8ff",
        3: "#91a7ff",
        4: "#748ffc",
        5: "#5c7cfa",
        6: "#4c6ef5",
        7: "#4263eb",
        8: "#3b5bdb",
        9: "#364fc7",
    },
    "blue": {
        0: "#e7f5ff",
        1: "#d0ebff",
        2: "#a5d8ff",
        3: "#74c0fc",
        4: "#4dabf7",
        5: "#339af0",
        6: "#228be6",
        7: "#1c7ed6",
        8: "#1971c2",
        9: "#1864ab",
    },
    "cyan": {
        0: "#e3fafc",
        1: "#c5f6fa",
        2: "#99e9f2",
        3: "#66d9e8",
        4: "#3bc9db",
        5: "#22b8cf",
        6: "#15aabf",
        7: "#1098ad",
        8: "#0c8599",
        9: "#0b7285",
    },
    "teal": {
        0: "#e6fcf5",
        1: "#c3fae8",
        2: "#96f2d7",
        3: "#63e6be",
        4: "#38d9a9",
        5: "#20c997",
        6: "#12b886",
        7: "#0ca678",
        8: "#099268",
        9: "#087f5b",
    },
    "green": {
        0: "#ebfbee",
        1: "#d3f9d8",
        2: "#b2f2bb",
        3: "#8ce99a",
        4: "#69db7c",
        5: "#51cf66",
        6: "#40c057",
        7: "#37b24d",
        8: "#2f9e44",
        9: "#2b8a3e",
    },
    "lime": {
        0: "#f4fce3",
        1: "#e9fac8",
        2: "#d8f5a2",
        3: "#c0eb75",
        4: "#a9e34b",
        5: "#94d82d",
        6: "#82c91e",
        7: "#74b816",
        8: "#66a80f",
        9: "#5c940d",
    },
    "yellow": {
        0: "#fff9db",
        1: "#fff3bf",
        2: "#ffec99",
        3: "#ffe066",
        4: "#ffd43b",
        5: "#fcc419",
        6: "#fab005",
        7: "#f59f00",
        8: "#f08c00",
        9: "#e67700",
    },
    "orange": {
        0: "#fff4e6",
        1: "#ffe8cc",
        2: "#ffd8a8",
        3: "#ffc078",
        4: "#ffa94d",
        5: "#ff922b",
        6: "#fd7e14",
        7: "#f76707",
        8: "#e8590c",
        9: "#d9480f",
    },
}

def main():
    event_number = 0
    nevents = 5
    e_cut = 0
    inFile = "/usatlas/u/atishelma/FCC/FCCAnalyses/2DPlotting/stage2/output.root"
    d = ROOT.RDataFrame("events", inFile)
     # can put some of the cuts we use to look at signal-like events
#     d = d.Filter("Zcand_m > 81") #sel1
    
    # now let's select the events
    d1 = d.Range(event_number, event_number + nevents) # Range([first, last[)
    d1 = d1.Define("H_phi", f"Hpart_phi[Hpart_energy>{e_cut}]")
    d1 = d1.Define("H_theta", f"Hpart_theta[Hpart_energy>{e_cut}]")
    d1 = d1.Define("H_energy", f"Hpart_energy[Hpart_energy>{e_cut}]")
    d1 = d1.Define("Z_phi", f"Zpart_phi[Zpart_energy>{e_cut}]")
    d1 = d1.Define("Z_theta", f"Zpart_theta[Zpart_energy>{e_cut}]")
    d1 = d1.Define("Z_energy", f"Zpart_energy[Zpart_energy>{e_cut}]")

    d1 = d1.Define("Jetconstituents_4_theta_1", f"jetconstituents_4_theta_1[jetconstituents_4_energy_1>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_theta_2", f"jetconstituents_4_theta_2[jetconstituents_4_energy_2>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_theta_3", f"jetconstituents_4_theta_3[jetconstituents_4_energy_3>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_theta_4", f"jetconstituents_4_theta_4[jetconstituents_4_energy_4>{e_cut}]")

    d1 = d1.Define("Jetconstituents_4_phi_1", f"jetconstituents_4_phi_1[jetconstituents_4_energy_1>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_phi_2", f"jetconstituents_4_phi_2[jetconstituents_4_energy_2>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_phi_3", f"jetconstituents_4_phi_3[jetconstituents_4_energy_3>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_phi_4", f"jetconstituents_4_phi_4[jetconstituents_4_energy_4>{e_cut}]")

    d1 = d1.Define("Jetconstituents_4_energy_1", f"jetconstituents_4_energy_1[jetconstituents_4_energy_1>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_energy_2", f"jetconstituents_4_energy_2[jetconstituents_4_energy_2>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_energy_3", f"jetconstituents_4_energy_3[jetconstituents_4_energy_3>{e_cut}]")
    d1 = d1.Define("Jetconstituents_4_energy_4", f"jetconstituents_4_energy_4[jetconstituents_4_energy_4>{e_cut}]")

    arr = d1.AsNumpy(["H_phi", "H_theta", "H_energy", "Z_phi", "Z_theta", "Z_energy", 
                      "truth_H_theta", "truth_H_phi", 
                      "truth_Zq1_theta", "truth_Zq1_phi",
                      "truth_Zq2_theta", "truth_Zq2_phi",
                      "truth_Hq1_theta", "truth_Hq1_phi",
                      "truth_Hq2_theta", "truth_Hq2_phi",
        "Jetconstituents_4_phi_1", "Jetconstituents_4_theta_1", "Jetconstituents_4_energy_1",
        "Jetconstituents_4_phi_2", "Jetconstituents_4_theta_2", "Jetconstituents_4_energy_2",
        "Jetconstituents_4_phi_3", "Jetconstituents_4_theta_3", "Jetconstituents_4_energy_3",
        "Jetconstituents_4_phi_4", "Jetconstituents_4_theta_4", "Jetconstituents_4_energy_4",
        "firstZ_m", "secondZ_m", "higgs_4_m",
        "firstZ_firstjet", "firstZ_secondjet", "secondZ_firstjet", "secondZ_secondjet",
        ])

    for idx in range (nevents):

        H_phi = np.asarray(arr['H_phi'][idx])
        H_theta = np.asarray(arr['H_theta'][idx])
        H_energy = np.asarray(arr['H_energy'][idx])
        Z_phi = np.asarray(arr['Z_phi'][idx])
        Z_theta = np.asarray(arr['Z_theta'][idx])
        Z_energy = np.asarray(arr['Z_energy'][idx])
        
        truth_H_theta = np.asarray(arr['truth_H_theta'][idx])
        truth_H_phi = np.asarray(arr['truth_H_phi'][idx])
        
        truth_Zq1_theta = np.asarray(arr['truth_Zq1_theta'][idx])
        truth_Zq1_phi = np.asarray(arr['truth_Zq1_phi'][idx])
        truth_Zq2_theta = np.asarray(arr['truth_Zq2_theta'][idx])
        truth_Zq2_phi = np.asarray(arr['truth_Zq2_phi'][idx])
        
        truth_Hq1_theta = np.asarray(arr['truth_Hq1_theta'][idx])
        truth_Hq1_phi = np.asarray(arr['truth_Hq1_phi'][idx])
        truth_Hq2_theta = np.asarray(arr['truth_Hq2_theta'][idx])
        truth_Hq2_phi = np.asarray(arr['truth_Hq2_phi'][idx])
        

        jc_phi_tmp = [np.asarray(arr['Jetconstituents_4_phi_1'][idx]),
                np.asarray(arr['Jetconstituents_4_phi_2'][idx]),
                np.asarray(arr['Jetconstituents_4_phi_3'][idx]),
                np.asarray(arr['Jetconstituents_4_phi_4'][idx])]
        jc_theta_tmp = [np.asarray(arr['Jetconstituents_4_theta_1'][idx]),
                np.asarray(arr['Jetconstituents_4_theta_2'][idx]),
                np.asarray(arr['Jetconstituents_4_theta_3'][idx]),
                np.asarray(arr['Jetconstituents_4_theta_4'][idx])]
        jc_energy_tmp = [np.asarray(arr['Jetconstituents_4_energy_1'][idx]),
                np.asarray(arr['Jetconstituents_4_energy_2'][idx]),
                np.asarray(arr['Jetconstituents_4_energy_3'][idx]),
                np.asarray(arr['Jetconstituents_4_energy_4'][idx])]

        mZ1 = arr['firstZ_m'][idx]
        mZ2 = arr['secondZ_m'][idx]
        higgsrecomass = arr['higgs_4_m'][idx]
        indices = [arr['firstZ_firstjet'][idx]-1, arr['firstZ_secondjet'][idx]-1,
                arr['secondZ_firstjet'][idx]-1, arr['secondZ_secondjet'][idx]-1]
        nconst = [len(jc) for jc in jc_phi_tmp]

        print("truth_Zq1_theta:",truth_Zq1_theta)
        print("truth_Zq1_phi:",truth_Zq1_phi)
        print("truth_Zq2_theta:",truth_Zq2_theta)
        print("truth_Zq2_phi:",truth_Zq2_phi)
        
        print("truth_Hq1_theta:",truth_Hq1_theta)
        print("truth_Hq1_phi:",truth_Hq1_phi)
        print("truth_Hq2_theta:",truth_Hq2_theta)
        print("truth_Hq2_phi:",truth_Hq2_phi)
        
        # trick to show correctly markers that overlap around pi in phi
        H_phi = np.concatenate((H_phi, H_phi+2*np.pi, H_phi-2*np.pi))
        H_theta = np.concatenate((H_theta, H_theta, H_theta))
        H_energy = np.concatenate((H_energy, H_energy, H_energy))
        Z_phi = np.concatenate((Z_phi, Z_phi+2*np.pi, Z_phi-2*np.pi))
        Z_theta = np.concatenate((Z_theta, Z_theta, Z_theta))
        Z_energy = np.concatenate((Z_energy, Z_energy, Z_energy))
        
        truth_H_phi = np.concatenate((truth_H_phi, truth_H_phi+2*np.pi, truth_H_phi-2*np.pi))
        truth_H_theta = np.concatenate((truth_H_theta, truth_H_theta, truth_H_theta))
        
        truth_Zq1_theta = np.concatenate((truth_Zq1_theta, truth_Zq1_theta, truth_Zq1_theta))
        truth_Zq2_theta = np.concatenate((truth_Zq2_theta, truth_Zq2_theta, truth_Zq2_theta))
        truth_Hq1_theta = np.concatenate((truth_Hq1_theta, truth_Hq1_theta, truth_Hq1_theta))
        truth_Hq2_theta = np.concatenate((truth_Hq2_theta, truth_Hq2_theta, truth_Hq2_theta))

        truth_Zq1_phi = np.concatenate((truth_Zq1_phi, truth_Zq1_phi+2*np.pi, truth_Zq1_phi-2*np.pi))
        truth_Zq2_phi = np.concatenate((truth_Zq2_phi, truth_Zq2_phi+2*np.pi, truth_Zq2_phi-2*np.pi))
        truth_Hq1_phi = np.concatenate((truth_Hq1_phi, truth_Hq1_phi+2*np.pi, truth_Hq1_phi-2*np.pi))
        truth_Hq2_phi = np.concatenate((truth_Hq2_phi, truth_Hq2_phi+2*np.pi, truth_Hq2_phi-2*np.pi))
        
        jc_phi = [np.concatenate((jc, jc+2*np.pi, jc-2*np.pi)) for jc in jc_phi_tmp]
        jc_theta = [np.concatenate((jc, jc, jc)) for jc in jc_theta_tmp]
        jc_energy = [np.concatenate((jc, jc, jc)) for jc in jc_energy_tmp]

        p = plt.rcParams
        #p["font.sans-serif"] = ["Nimbus Sans"]
#         p["font.sans-serif"] = ["Liberation Sans"]
        p["mathtext.fontset"] = "stixsans"
        #p["font.weight"] = "light"

        fig = plt.figure(figsize=(8,5), tight_layout=True, dpi=200)
        ax = plt.subplot(1,1,1, aspect=1, xlim = [-np.pi, np.pi], ylim = [0, 3.8])

        ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
        ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "$0$", r"$+\pi/2$", r"$+\pi$"])
        ax.set_yticks([0, np.pi/2, np.pi])
        ax.set_yticklabels(["$0$", r"$+\pi/2$", r"$+\pi$"])
        ax.set_xlabel(r"$\mathbf{\phi}$")
        ax.set_ylabel(r"$\mathbf{\theta}$")
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.spines.left.set_position(("data", -np.pi-0.2))
        ax.spines.left.set_bounds(0, np.pi)
        ax.spines.bottom.set_position(("data", -0.2))

        # use red-ish colors for jets used to create on-shell Z
        # use green-ish colors for jets used to create off-shell Z
        jet_colors = [colors['red'][4], colors['grape'][4],
                colors['lime'][4], colors['teal'][4]]

        js = []

        for i, ii in enumerate(indices):
            js.append(ax.scatter(jc_phi[ii], jc_theta[ii], s=100*jc_energy[ii]**.5,
                c=jet_colors[i], alpha=0.4,
                label=f"Jet {ii+1}, {nconst[ii]} constituents"))

        print("len(H_phi):",len(H_phi))
        print("len(H_theta):",len(H_theta))
        
        print("len(Z_phi):",len(Z_phi))
        print("len(Z_theta):",len(Z_theta))
        
        print("len(truth_H_phi):",len(truth_H_phi))
        print("len(truth_H_theta):",len(truth_H_theta))

        H_truth=ax.scatter(H_phi, H_theta, s=10*H_energy**.5,
                c=colors['orange'][8], alpha=1, marker='P',
                label="Truth particles from Higgs(bb)")
        Z_truth=ax.scatter(Z_phi, Z_theta, s=10*Z_energy**.5,
                c=colors['green'][8], alpha=1, marker='P',
                label="Truth particles from Z(cc)")

        ax.legend(handles=[H_truth, js[0], js[1], Z_truth, js[2], js[3]], loc='upper left', ncol=2)
#         ax.text(0, 1.0, fr"Reconstructed masses: $m(Z_1) = {mZ1:.1f}$ GeV, $m(Z_2) = {mZ2:.1f}$ GeV, $m(H) = {higgsrecomass:.1f}$ GeV",
        ax.text(0, 1.0, fr"",
            transform=ax.transAxes)
        ax.set_title("Z(cc)H(bb) events, after selections, Durham-kt N=4", loc='left',
                fontdict={'fontweight':'bold'}, y=1.03)

        plt.savefig(f"outputs/Display_{idx}.png", dpi=500)
        plt.close()
    pass


if __name__ == "__main__":
    main()
