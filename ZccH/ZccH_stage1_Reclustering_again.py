"""
25 January 2023
Abraham Tishelman-Charny

The purpose of this python module is to run stage 1 of the FCC analysis for Z(cc)H. Started with examples in repo. 

"""

import os 
import urllib.request

# ____________________________________________________________
def get_file_path(url, filename):
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        urllib.request.urlretrieve(url, os.path.basename(url))
        return os.path.basename(url)

# ____________________________________________________________

processList = {
    # Z(cc)H by higgs final state 
    # 'wzp6_ee_ccH_HWW_ecm240':{'chunks':20},
    # 'wzp6_ee_ccH_Hgg_ecm240' : {'chunks':20},
    # 'wzp6_ee_ccH_HZa_ecm240' : {'chunks':20},
    # 'wzp6_ee_ccH_Hss_ecm240' : {'chunks':20},
    # 'wzp6_ee_ccH_Hcc_ecm240' : {'chunks':20},
    # 'wzp6_ee_ccH_Hmumu_ecm240':{'chunks':20},
    # 'wzp6_ee_ccH_HZZ_ecm240' : {'chunks':20},	
    # 'wzp6_ee_ccH_Htautau_ecm240' : {'chunks':20},
    # 'wzp6_ee_ccH_Haa_ecm240' : {'chunks':20},
    'wzp6_ee_ccH_Hbb_ecm240':{'chunks':20},

    # backgrounds
    #'p8_ee_WW_ecm240' : {'chunks':20},
    #'p8_ee_ZZ_ecm240' : {'chunks':20}
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/" 

## input file needed for unit test in CI
testFile = "https://fccsw.web.cern.ch/fccsw/testsamples/wzp6_ee_nunuH_Hss_ecm240.root"

## latest particle transformer model, trainied on 9M jets in winter2023 samples
model_name = "fccee_flavtagging_edm4hep_wc_v1"

## model files needed for unit testing in CI
url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)

## model files locally stored on /eos
model_dir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

## get local file, else download from url
weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)


from addons.ONNXRuntime.python.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.python.jetClusteringHelper import ExclusiveJetClusteringHelper

jetFlavourHelper = None
jetClusteringHelper = None

# Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis:
    # __________________________________________________________
    # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        global jetClusteringHelper
        global jetFlavourHelper

        from examples.FCCee.weaver.config import collections, njets

        tag = ""

        ## define jet clustering parameters
        jetClusteringHelper = ExclusiveJetClusteringHelper(collections["PFParticles"], njets, tag)

        ## run jet clustering
        df = jetClusteringHelper.define(df)

        ## define jet flavour tagging parameters

        jetFlavourHelper = JetFlavourHelper(
            collections,
            jetClusteringHelper.jets,
            jetClusteringHelper.constituents,
            tag,
        )

        ## define observables for tagger
        df = jetFlavourHelper.define(df)

        ## tagger inference
        df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df)

        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():

        ##  outputs jet properties
        branchList = jetClusteringHelper.outputBranches()

        ## outputs jet scores and constituent breakdown
        branchList += jetFlavourHelper.outputBranches()

        return branchList


# #Mandatory: RDFanalysis class where the use defines the operations on the TTree
# class RDFanalysis():

#     #__________________________________________________________
#     #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
#     def analysers(df):
#         global jetClusteringHelper
#         global jetFlavourHelper
        
        
        # df2 = (
        #     df
        #        .Alias("Jet3","Jet#3.index")
        #        .Define("DeltaR", "2.0")
        #        #define the RP px, py, pz and e
        #        .Define("RP_px",          "ReconstructedParticle::get_px(ReconstructedParticles)")
        #        .Define("RP_py",          "ReconstructedParticle::get_py(ReconstructedParticles)")
        #        .Define("RP_pz",          "ReconstructedParticle::get_pz(ReconstructedParticles)")
        #        .Define("RP_e",           "ReconstructedParticle::get_e(ReconstructedParticles)")
        #        .Define("RP_m",           "ReconstructedParticle::get_mass(ReconstructedParticles)")
        #        .Define("RP_q",           "ReconstructedParticle::get_charge(ReconstructedParticles)")

        #        #build pseudo jets with the RP, using the interface that takes px,py,pz,m for better
        #        #handling of rounding errors
        #        .Define("pseudo_jets",    "JetClusteringUtils::set_pseudoJets_xyzm(RP_px, RP_py, RP_pz, RP_m)")
        #        #.Define("pseudo_jets",    "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

        #        #run jet clustering with all reconstructed particles. kt_algorithm, R=4, exclusive clustering, exactly 4 jets, E0-scheme
        #        .Define("FCCAnalysesJets_kt", "JetClustering::clustering_kt(DeltaR, 2, 4, 0, 10)(pseudo_jets)") # first param is deltaR
        #        #get the jets out of the struct
        #        .Define("jets_kt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_kt)")
        #        #get the jets constituents out of the struct
        #        .Define("jetconstituents_kt","JetClusteringUtils::get_constituents(FCCAnalysesJets_kt)")
        #        #get some variables
        #        .Define("jets_kt_e",        "JetClusteringUtils::get_e(jets_kt)")
        #        .Define("jets_kt_px",        "JetClusteringUtils::get_px(jets_kt)")
        #        .Define("jets_kt_py",        "JetClusteringUtils::get_py(jets_kt)")
        #        .Define("jets_kt_pz",        "JetClusteringUtils::get_pz(jets_kt)")
        #        .Define("jets_kt_m",        "JetClusteringUtils::get_m(jets_kt)")

        #        #run jet clustering with all reconstructed particles. ee_genkt_algorithm, R=0.5, inclusive clustering, E-scheme
        #        .Define("FCCAnalysesJets_ee_genkt", "JetClustering::clustering_ee_genkt(DeltaR, 0, 0, 0, 0, -1)(pseudo_jets)")
        #        #get the jets out of the struct
        #        .Define("jets_ee_genkt",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt)")
        #        #get the jets constituents out of the struct
        #        .Define("jetconstituents_ee_genkt","JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt)")
        #        #get some variables
        #        .Define("jets_ee_genkt_px",        "JetClusteringUtils::get_px(jets_ee_genkt)")
        #        .Define("jets_ee_genkt_py",        "JetClusteringUtils::get_py(jets_ee_genkt)")
        #        .Define("jets_ee_genkt_pz",        "JetClusteringUtils::get_pz(jets_ee_genkt)")
        #        .Define("jets_ee_genkt_m",        "JetClusteringUtils::get_m(jets_ee_genkt)")

        #        #run jet clustering with all reconstructed particles. valencia_algorithm, R=0.5, inclusive clustering, E-scheme
        #        .Define("FCCAnalysesJets_valencia", "JetClustering::clustering_valencia(DeltaR, 1, 6, 0, 0, 1., 1.)(pseudo_jets)")

        #        #get the jets out of the struct
        #        .Define("jets_valencia",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_valencia)")
        #        #get the jets constituents out of the struct
        #        .Define("jetconstituents_valencia","JetClusteringUtils::get_constituents(FCCAnalysesJets_valencia)")
        #        #get some variables
        #        .Define("jets_valencia_px",        "JetClusteringUtils::get_px(jets_valencia)")
        #        .Define("jets_valencia_py",        "JetClusteringUtils::get_py(jets_valencia)")
        #        .Define("jets_valencia_pz",        "JetClusteringUtils::get_pz(jets_valencia)")
        #        .Define("jets_valencia_m",        "JetClusteringUtils::get_m(jets_valencia)")

        #        #run jet clustering with all reconstructed particles. jade_algorithm, R=0.5, exclusive clustering, exactly 4 jets, sorted by E, E0-scheme
        #        .Define("FCCAnalysesJets_jade", "JetClustering::clustering_jade(DeltaR, 2, 4, 1, 10)(pseudo_jets)")

        #        #get the jets out of the struct
        #        .Define("jets_jade",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_jade)")
        #        #get the jets constituents out of the struct
        #        .Define("jetconstituents_jade","JetClusteringUtils::get_constituents(FCCAnalysesJets_jade)")
        #        #get some variables
        #        .Define("jets_jade_px",        "JetClusteringUtils::get_px(jets_jade)")
        #        .Define("jets_jade_py",        "JetClusteringUtils::get_py(jets_jade)")
        #        .Define("jets_jade_pz",        "JetClusteringUtils::get_pz(jets_jade)")
        #        .Define("jets_jade_m",        "JetClusteringUtils::get_m(jets_jade)")
        #        .Define("jets_jade_flavour",   "JetTaggingUtils::get_flavour(jets_jade, Particle)")
        #        .Define("jets_jade_btag",      "JetTaggingUtils::get_btag(jets_jade_flavour, 0.80)")
        #        .Define("jets_jade_btag_true", "JetTaggingUtils::get_btag(jets_jade_flavour, 1.0)")
        #        .Define("jets_jade_ctag",      "JetTaggingUtils::get_ctag(jets_jade_flavour, 0.10)")
        #        .Define("jets_jade_ctag_true",      "JetTaggingUtils::get_ctag(jets_jade_flavour, 1.0)")

        #        .Define("JET_btag",       "ReconstructedParticle::getJet_btag(Jet3, ParticleIDs, ParticleIDs_0)")
        #        .Define("EVT_nbtag",      "ReconstructedParticle::getJet_ntags(JET_btag)")

        #        .Define('EVT_thrust',     'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
        #        .Define('RP_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')
        #        .Define('EVT_thrust_x',   "EVT_thrust.at(0)")
        #        .Define('EVT_thrust_y',   "EVT_thrust.at(1)")
        #        .Define('EVT_thrust_z',   "EVT_thrust.at(2)")
        #        .Define('EVT_thrust_val', "EVT_thrust.at(3)")

        #        .Define('EVT_sphericity',     'Algorithms::minimize_sphericity("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
        #        .Define('EVT_sphericity_x',   "EVT_sphericity.at(0)")
        #        .Define('EVT_sphericity_y',   "EVT_sphericity.at(1)")
        #        .Define('EVT_sphericity_z',   "EVT_sphericity.at(2)")
        #        .Define('EVT_sphericity_val', "EVT_sphericity.at(3)")
        #        .Define('RP_sphericityangle', 'Algorithms::getAxisCosTheta(EVT_sphericity, RP_px, RP_py, RP_pz)')

        #        .Define('RP_hemis0_mass',   "Algorithms::getAxisMass(0)(RP_thrustangle, RP_e, RP_px, RP_py, RP_pz)")
        #        .Define('RP_hemis1_mass',   "Algorithms::getAxisMass(1)(RP_thrustangle, RP_e, RP_px, RP_py, RP_pz)")

        #        .Define("RP_total_mass",    "Algorithms::getMass(ReconstructedParticles)")

        # )
        # return df2
    






    # #__________________________________________________________
    # #Mandatory: output function, please make sure you return the branchlist as a python list
    # def output():
    #     branchList = [
    #             "RP_px", "RP_py", "RP_pz", "RP_e", "RP_m", "RP_q",

    #             "JET_btag",
    #             "EVT_nbtag",

    #             "EVT_thrust_x", "EVT_thrust_y", "EVT_thrust_z", "EVT_thrust_val",

    #             "EVT_sphericity_x", "EVT_sphericity_y", "EVT_sphericity_z", "EVT_sphericity_val",

    #             "RP_thrustangle",
    #             "RP_sphericityangle",

    #             "RP_hemis0_mass",
    #             "RP_hemis1_mass",
    #             "RP_total_mass",

    #             "jets_kt_e",
    #             "jets_kt_px",
    #             "jets_kt_py",
    #             "jets_kt_pz",
    #             "jets_kt_m",
    #             "jetconstituents_kt",

    #             "jets_ee_genkt_px",
    #             "jets_ee_genkt_py",
    #             "jets_ee_genkt_pz",
    #             "jets_ee_genkt_m",
    #             "jetconstituents_ee_genkt",

    #             "jets_valencia_px",
    #             "jets_valencia_py",
    #             "jets_valencia_pz",
    #             "jets_valencia_m",
    #             "jetconstituents_valencia",

    #             "jets_jade_px",
    #             "jets_jade_py",
    #             "jets_jade_pz",
    #             "jets_jade_m",
    #             "jets_jade_ctag",
    #             "jets_jade_ctag_true",
    #             "jets_jade_btag",
    #             "jets_jade_btag_true",
    #             "jetconstituents_jade",

    #             ]
    #     return branchList