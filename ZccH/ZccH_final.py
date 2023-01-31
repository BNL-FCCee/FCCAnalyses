"""
25 January 2023
Abraham Tishelman-Charny

The purpose of this python module is to run the 'final' step of the FCC analysis for Z(cc)H. Started with examples in repo. 

Does batch mode not work for the final step?

"""

import yaml 

configFile = "/afs/cern.ch/work/a/atishelm/private/FCCAnalyses/ZccH/RunConfig.yaml"
with open(configFile, 'r') as cfg:
    values = yaml.safe_load(cfg)
    
    batch = values["batch"]
    EOSoutput = values["EOSoutput"]
    JobName = values["JobName"]

print("batch:",batch)
print("EOSoutput:",EOSoutput)
print("JobName:",JobName)

processList = {
    # Z(cc)H by higgs final state 
    'wzp6_ee_ccH_HWW_ecm240':{'chunks':20},
    'wzp6_ee_ccH_Hgg_ecm240' : {'chunks':20},
    'wzp6_ee_ccH_HZa_ecm240' : {'chunks':20},
    'wzp6_ee_ccH_Hss_ecm240' : {'chunks':20},
    'wzp6_ee_ccH_Hcc_ecm240' : {'chunks':20},
    'wzp6_ee_ccH_Hmumu_ecm240':{'chunks':20},
    'wzp6_ee_ccH_HZZ_ecm240' : {'chunks':20},	
    'wzp6_ee_ccH_Htautau_ecm240' : {'chunks':20},
    'wzp6_ee_ccH_Haa_ecm240' : {'chunks':20},
    'wzp6_ee_ccH_Hbb_ecm240':{'chunks':20},

    # backgrounds
    #'p8_ee_WW_ecm240' : {'chunks':20},
    #'p8_ee_ZZ_ecm240' : {'chunks':20}
}

#Link to the dictonary that contains all the cross section informations etc...
procDict = "FCCee_procDict_winter2023_IDEA.json" # do we have one for winter2023?

if(EOSoutput):
    # inputDir    = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/stage2/"
    inputDir    = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/stage1/"
    outputDir = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/final/"
    eosType = "eosuser"
else:
    inputDir    = f"{JobName}/stage2/"
    outputDir   = f"{JobName}/final/"

nCPUS       = 4
runBatch    = batch
batchQueue = "longlunch"
#compGroup = "group_u_FCC.local_gen"

#Add MySample_p8_ee_ZH_ecm240 as it is not an offical process
#procDictAdd={"MySample_p8_ee_ZH_ecm240":{"numberOfEvents": 10000000, "sumOfWeights": 10000000, "crossSection": 0.201868, "kfactor": 1.0, "matchingEfficiency": 1.0}}

#Number of CPUs to use
nCPUS = 4

#produces ROOT TTrees, default is False
doTree = True

###Dictionnay of the list of cuts. The key is the name of the selection that will be added to the output file
# cutList = {"sel0":"Zcand_q == 0",
#             "sel1":"Zcand_q == -1 || Zcand_q == 1",
#             "sel2":"Zcand_m > 80 && Zcand_m < 100",
#             "sel3":"MyFilter==true && (Zcand_m < 80 || Zcand_m > 100)"
#             }

cutList = {
    "sel0" : "1",
    # "sel1" : "selected_jets_pt > 20",
    # "sel2" : "selected_jets_pt > 30",
    # "sel3" : "selected_jets_pt > 40",
    # "sel4" : "selected_jets_pt > 50",

    #"sel0" : "Zcand_m > 40 && Zcand_m < 120",
}


#Dictionary for the ouput variable/hitograms. The key is the name of the variable in the output files. "name" is the name of the variable in the input file, "title" is the x-axis label of the histogram, "bin" the number of bins of the histogram, "xmin" the minimum x-axis value and "xmax" the maximum x-axis value.
histoList = {

            # # jets 
            # "N_selected_jets",
            # "selected_jets_pt",
            # "selected_jets_y",
            # "selected_jets_p",
            # "selected_jets_e",

            # # dijets pairs 
            # "N_dijet_pair",
            # "dijet_pair_all_pt",
            # "dijet_pair_all_m",
            # "dijet_pair_all_recoil_m",

            # # Zed candidate near z mass peak 
            # "dijet_pair_nearZpeak_m",           
            # "dijet_pair_nearZpeak_pt",           
            # "dijet_pair_nearZpeak_recoil_m",     

    "N_selected_jets": {"name":"N_selected_jets", "title":"N_{jets}", "bin":10, "xmin": 0, "xmax" : 10},
    "jets_pt" : {"name":"selected_jets_pt","title":"All Jets p_{T} [GeV]","bin":125,"xmin":0,"xmax":300},
    "jets_y" : {"name":"selected_jets_y","title":"All Jets y [rad]","bin":40,"xmin":-10,"xmax":10},
    "jets_p" : {"name":"selected_jets_p","title":"All Jets p [GeV]","bin":125,"xmin":0,"xmax":300},
    "jets_e" : {"name":"selected_jets_e","title":"All Jets e [GeV]","bin":125,"xmin":0,"xmax":300},

    "N_dijet_pair": {"name":"N_dijet_pair", "title":"N_{dijets}", "bin":10, "xmin": 0, "xmax" : 10},
    "dijet_pair_all_pt" : {"name":"dijet_pair_all_pt","title":"All Jets p_{T} [GeV]","bin":125,"xmin":0,"xmax":300},
    "dijet_pair_all_m" : {"name":"dijet_pair_all_m","title":"m_{jj} [GeV]","bin":125,"xmin":0,"xmax":250},
    "dijet_pair_all_recoil_m" : {"name":"dijet_pair_all_recoil_m","title":"recoil mass [GeV]","bin":125,"xmin":0,"xmax":250},

    "dijet_pair_nearZpeak_pt" : {"name":"dijet_pair_nearZpeak_pt","title":"dijet_pair_nearZpeak_pt [GeV]","bin":125,"xmin":0,"xmax":300},
    "dijet_pair_nearZpeak_m" : {"name":"dijet_pair_nearZpeak_m","title":"m_{jj} (nearZpeak) [GeV]","bin":125,"xmin":0,"xmax":250},
    "dijet_pair_nearZpeak_recoil_m" : {"name":"dijet_pair_nearZpeak_recoil_m","title":"recoil mass (nearZpeak) [GeV]","bin":125,"xmin":0,"xmax":250},

    # "jets_pt_0" : {"name":"selected_jets_pt_0","title":"Jet 0 p_{T} [GeV]","bin":125,"xmin":0,"xmax":300},
    # "jets_y_0" : {"name":"selected_jets_y_0","title":"Jet 0 y [rad]","bin":40,"xmin":-10,"xmax":10},
    # "jets_p_0" : {"name":"selected_jets_p_0","title":"Jet 0 p [GeV]","bin":125,"xmin":0,"xmax":300},
    # "jets_e_0" : {"name":"selected_jets_e_0","title":"Jet 0 e [GeV]","bin":125,"xmin":0,"xmax":300},    
    
    # "mjj":{"name":"Zcand_m","title":"m_{jj} [GeV]","bin":125,"xmin":0,"xmax":250},
    # "hadronic_recoil_m":{"name":"Zcand_recoil_m","title":"Z hadronic recoil [GeV]","bin":100,"xmin":0,"xmax":400},
}
