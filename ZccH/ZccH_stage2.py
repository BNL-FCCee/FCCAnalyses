"""
25 January 2023
Abraham Tishelman-Charny

The purpose of this python module is to run stage 2 of the FCC analysis for Z(cc)H. Started with examples in repo. 

"""

import yaml 

configFile = "/afs/cern.ch/work/a/atishelm/private/FCCAnalyses/ZccH/RunConfig.yaml"
with open(configFile, 'r') as cfg:
    values = yaml.safe_load(cfg)
    
    batch = values["batch"]
    EOSoutput = values["EOSoutput"]

print("batch:",batch)
print("EOSoutput:",EOSoutput)

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
    'p8_ee_WW_ecm240' : {'chunks':20},
    'p8_ee_ZZ_ecm240' : {'chunks':20}
}

if(EOSoutput):
    inputDir    = "/eos/user/a/atishelm/ntuples/FCC/ZccH/stage1/"
    outputDir = "/eos/user/a/atishelm/ntuples/FCC/ZccH/stage2/"
    eosType = "eosuser"
else:
    inputDir    = "ZccH/stage1/"
    outputDir   = "ZccH/stage2/"

nCPUS       = 4
runBatch    = batch
batchQueue = "longlunch"
#compGroup = "group_u_FCC.local_gen"


#USER DEFINED CODE
# import ROOT
# ROOT.gInterpreter.Declare("""
# bool myFilter(ROOT::VecOps::RVec<float> mass) {
#     for (size_t i = 0; i < mass.size(); ++i) {
#         if (mass.at(i)>1. && mass.at(i)<300.)
#             return true;
#     }
#     return false;
# }
# """)
#END USER DEFINED CODE

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (df
               #Filter to have exactly one Zcc candidate
               .Filter("zed_hadronic_m.size() == 1") # exactly one 

               # create branches with jet (transverse momentum, rapidity, total momentum, energy)
               .Define("selected_jets_pt_0", "selected_jets_pt[0]") # what does [0] take? leading? 
               .Define("selected_jets_y_0", "selected_jets_y[0]")
               .Define("selected_jets_p_0", "selected_jets_p[0]")
               .Define("selected_jets_e_0", "selected_jets_e[0]")

               #Define Z candidate mass
               .Define("Zcand_m","zed_hadronic_m[0]")
               #Define Z candidate recoil mass
               .Define("Zcand_recoil_m","zed_hadronic_recoil_m[0]")
               #Define Z candidate pt
               .Define("Zcand_pt","zed_hadronic_pt[0]")
               #Define new var rdf entry (example)
               #.Define("entry", "rdfentry_")
               #Define a weight based on entry (inline example of possible operations)
               #.Define("weight", "return 1./(entry+1)")
               #Define a variable based on a custom filter
               #.Define("MyFilter", "myFilter(zed_leptonic_m)") # maybe choose flavour here. 
               )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list.
    def output():
        branchList = [
            "selected_jets_pt",
            "selected_jets_y",
            "selected_jets_p",
            "selected_jets_e",            
            "selected_jets_pt_0",
            "selected_jets_y_0",
            "selected_jets_p_0",
            "selected_jets_e_0",            
            "Zcand_m", 
            "Zcand_pt", 
            "Zcand_recoil_m",
            #"entry",
            #"weight",
        ]
        return branchList
