"""
25 January 2023
Abraham Tishelman-Charny

The purpose of this python module is to run stage 1 of the FCC analysis for Z(cc)H. Started with examples in repo. 

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

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/" 

if(EOSoutput):
    outputDir = "/eos/user/a/atishelm/ntuples/FCC/ZccH/stage1/"
    eosType = "eosuser"
else:
    outputDir   = "ZccH/stage1/"

nCPUS       = 4
runBatch    = batch
batchQueue = "longlunch"
#compGroup = "group_u_FCC.local_gen"

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df
            # define an alias for jet index collection
            .Alias("Jet2", "Jet#2.index")

            # define the jet collection 
            .Define("jets", "ReconstructedParticle::get(Jet2, ReconstructedParticles)")

            # define collections for jet flavour, and mother particle? 

            # select jets on pT 
            .Define("selected_jets", "ReconstructedParticle::sel_pt(10.)(jets)") # 10 GeV selection on jets 

            # create branches with jet (transverse momentum, rapidity, total momentum, energy)
            .Define("selected_jets_pt", "ReconstructedParticle::get_pt(selected_jets)")
            .Define("selected_jets_y", "ReconstructedParticle::get_y(selected_jets)")
            .Define("selected_jets_p", "ReconstructedParticle::get_p(selected_jets)")
            .Define("selected_jets_e", "ReconstructedParticle::get_e(selected_jets)")

            # find zed candidates from di-jet resonances
            .Define("zed_hadronic", "ReconstructedParticle::resonanceBuilder(91)(selected_jets)")

            # create branches with zed (mass, transverse momentum, recoil, recoil mass, charge)
            .Define("zed_hadronic_m", "ReconstructedParticle::get_mass(zed_hadronic)")
            .Define("zed_hadronic_pt", "ReconstructedParticle::get_pt(zed_hadronic)")
            .Define("zed_hadronic_recoil", "ReconstructedParticle::recoilBuilder(240)(zed_hadronic)")
            .Define("zed_hadronic_recoil_m", "ReconstructedParticle::get_mass(zed_hadronic_recoil)")

            # Filter at least one candidate 
            .Filter("zed_hadronic_recoil_m.size()>0")

        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            "selected_jets_pt",
            "selected_jets_y",
            "selected_jets_p",
            "selected_jets_e",
            "zed_hadronic_pt",
            "zed_hadronic_m",
            "zed_hadronic_recoil_m"
        ]
        return branchList
