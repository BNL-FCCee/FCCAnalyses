"""
January 2023
Abraham Tishelman-Charny
The purpose of this python module is to perform initial selections and variable definitions for processing FCC files.

July 2024
Anna Elizabeth Connelly

"""

import ROOT
import os
import urllib.request
#import yaml 
import sys


sys.path.append("/usatlas/u/aconnelly/IzaFCCAnalysis")

from examples.FCCee.weaver.config import collections
from CustomDefinitions import CustomDefinitions
from examples.FCCee.weaver.config import (
    #variables_pfcand,
    variables_jet
)
from addons.ONNXRuntime.python.jetFlavourHelper import JetFlavourHelper
#Exclusive Clustering Class
from addons.FastJet.python.jetClusteringHelper import ExclusiveJetClusteringHelper
#Inclusive Clustering Class
from addons.FastJet.python.jetClusteringHelper import InclusiveJetClusteringHelper

# originally had YAML config here. Not strictly necessary. Check previous commits if you want an example.
# batch = 1 # use HTCondor
# EOSoutput = 0 # output to EOS
# JobName = "ZHadronic_4JetReco" # job named used for output directory

#User input variables:
#number of jets in exclusive reclustering 
njets = 4 

#radius in inclusive clustering
rad = 1.0

#set algorithms-- inclusive algorithms -- 0-antikt, 1-inclusive eekt  2-cambridge
alg = 0 

#set sorted -- for inclusive-- 0-sorted by pt, 1-sorted by E
sort = 1

#energy cut to PseudoJets
ecut = 10

#variables used for reference in other files
vars = [njets, rad, alg, sort, ecut]

outputDir  = "/usatlas/u/aconnelly/IzaFCCAnalysis/zhanalysis/root"

processList = {

    
    }

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag = "FCCee/winter2023/IDEA/" 
procDict = "FCCee_procDict_winter2023_IDEA.json" 

# if(EOSoutput):
#     # example output directory on EOS - note that by default this includes copying the output file from one location to another
#     #outputDirEos = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/stage1/" # if you define outputDirEos, this process creates the file locally and copies it to eos.
#     #outputDir = f"/eos/user/a/atishelm/ntuples/FCC/ZH_Hadronic_4JetReco/"

#     eosType = "eosuser" # specify as necessary

#runBatch   = False
# batchQueue = "testmatch" 

# # Define any functionality which is not implemented in FCCAnalyses

# # ____________________________________________________________
def get_file_path(url, filename):
    print("Looking for file:",filename)
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        urllib.request.urlretrieve(url, os.path.basename(url))
        return os.path.basename(url)

# # ____________________________________________________________

ROOT.gInterpreter.Declare(CustomDefinitions)

## input file needed for unit test in CI
# testFile = "https://fccsw.web.cern.ch/fccsw/testsamples/wzp6_ee_nunuH_Hss_ecm240.root"

# ## latest particle transformer model, trained on 9M jets in winter2023 samples - need to separate train/test samples?
model_name = "fccee_flavtagging_edm4hep_wc"

# ## model files needed for unit testing in CI
url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
#url_model_dir = "/usatlas/u/aconnelly/IzaFCCAnalysis"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
#url_preproc = "/usatlas/u/aconnelly/IzaFCCAnalysis/fccee_flavtagging_edm4hep_wc.json"
url_model = "{}/{}.onnx".format(url_model_dir, model_name)
#url_model = "/usatlas/u/aconnelly/IzaFCCAnalysis/fccee_flavtagging_edm4hep_wc.json"

# ## model files locally stored on /eos
# if(batch):
#     #model_dir = "/usatlas/u/ivelisce/FCC_at_BNL/FCCAnalyses/"
#     model_dir = "/usatlas/u/aconnelly/IzaFCCAnalysis/"
# else: model_dir = "./"

model_dir = "/usatlas/u/aconnelly/IzaFCCAnalysis"
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

# ## get local file, else download from url
weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)

# ee_ktClustering = None
# # ee_ktFlavourHelper = None
antiktClustering = None
antiktFlavourHelper = None

def analysis_sequence(df):
    collections["Electrons"] = "Electron"
    collections["Muons"] = "Muon"
    collections["Photons"] = "Photons"

    df = (
        # electrons
        df.Alias("Electron0", "{}#0.index".format(collections["Electrons"]))
        .Define(
            "electrons",
            "ReconstructedParticle::get(Electron0, {})".format(collections["PFParticles"]),
        )
        .Define("event_nel", "electrons.size()")
        .Define("electrons_p", "ReconstructedParticle::get_p(electrons)[0]")
        
        # muons
        .Alias("Muon0", "{}#0.index".format(collections["Muons"]))
        .Define(
            "muons",
            "ReconstructedParticle::get(Muon0, {})".format(collections["PFParticles"]),
        )
        .Define("event_nmu", "muons.size()")
        .Define("muons_p", "ReconstructedParticle::get_p(muons)[0]")
        
        #Get kinematics variables needed for selection later
        .Define("P4_vis", "ReconstructedParticle::get_P4vis({})".format(collections["PFParticles"]))
        .Define("vis_M", "P4_vis.M()")
        .Define("vis_E", "P4_vis.E()")
        .Define("P3_vis","TVector3(P4_vis.Px(), P4_vis.Py(), P4_vis.Pz())")
        .Define("vis_theta", "P3_vis.Theta()")

        #EVENTWIDE VARIABLES: Access quantities that exist only once per event, such as the missing energy (despite the name, 
        # the MissingET collection contains the total missing energy)
        .Define("RecoMissingEnergy_e", "ReconstructedParticle::get_e(MissingET)")
        .Define("RecoMissingEnergy_p", "ReconstructedParticle::get_p(MissingET)")
        .Define("RecoMissingEnergy_pt", "ReconstructedParticle::get_pt(MissingET)")
        .Define("RecoMissingEnergy_px", "ReconstructedParticle::get_px(MissingET)") #x-component of RecoMissingEnergy
        .Define("RecoMissingEnergy_py", "ReconstructedParticle::get_py(MissingET)") #y-component of RecoMissingEnergy
        .Define("RecoMissingEnergy_pz", "ReconstructedParticle::get_pz(MissingET)") #z-component of RecoMissingEnergy
        .Define("RecoMissingEnergy_eta", "ReconstructedParticle::get_eta(MissingET)")
        .Define("RecoMissingEnergy_theta", "ReconstructedParticle::get_theta(MissingET)")
        .Define("RecoMissingEnergy_phi", "ReconstructedParticle::get_phi(MissingET)") #angle of RecoMissingEnergy
    )
    return df


#def jet_sequence(df, njets, exclusive):
#def jet_sequence(df, njets, rad, alg):

def jet_sequence(df,rad, alg, sort, ecut):

    # global ee_ktClustering
    # global ee_ktFlavourHelper
    global antiktClustering
    global antiktFlavourHelper

    ##First inclusive algorithm clustering --- Antikt
    tag = ""

    ## define jet clustering parameters
    antiktClustering = InclusiveJetClusteringHelper(collections["PFParticles"],rad, alg, sort, ecut, tag)
  
    ## runs inclusive antikt jet clustering 
    #extract all jet observables from pseudojets

    df = antiktClustering.define(df)

    ## define jet flavour tagging parameters
    antiktFlavourHelper = JetFlavourHelper(
        collections,
        antiktClustering.jets,
        antiktClustering.constituents,
        tag,
    )

    ## define observables for tagger
    #antiktClustering.jets refers to the pseudojets from get_pseudojets from the clustering algorithm
    #converting the pseudojets to TLorentz vectors 
    df = antiktFlavourHelper.define(df)
    df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets({})".format("antiktClustering.jets"))


    # apply energy correction
    jet_reco_vars = ["e", "p", "px", "py", "pz", "m", "theta"]

    for jet_reco_var in jet_reco_vars:
        df=(df.Define("recojet_{}".format(jet_reco_var), "JetClusteringUtils::get_{}(jet)".format(jet_reco_var)))
    
    # phi has slightly different naming
    df=(df.Define("recojet_phi", "JetClusteringUtils::get_phi_std(jet)"))
 
    ###
    df = df.Define("jets_tlv_corr", "FCCAnalyses::energyReconstructFourJet(recojet_px, recojet_py, recojet_pz, recojet_e)")

    jet_corr_vars = ["e", "px", "py", "pz"]
    for jet_corr_var in jet_corr_vars: 
         df = df.Define("jet_{}_corr".format(jet_corr_var), "FCCAnalyses::TLVHelpers::get_{}(jets_tlv_corr)".format(jet_corr_var))

    df = df.Define("all_invariant_masses", "JetConstituentsUtils::all_invariant_masses(jet_p4)")
    df = df.Define("recoil_masses", "all_recoil_masses(jet_p4)")
    
    ## tagger inference
    df = antiktFlavourHelper.inference(weaver_preproc, weaver_model, df) 

    ## define variables using tagger inference outputs
    df = df.Define("recojetpair_isC", "SumFlavorScores(recojet_isC)") 
    df = df.Define("recojetpair_isB", "SumFlavorScores(recojet_isB)") 

    df = df.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(_jet)")
    df = df.Define("jets_truth", "FCCAnalyses::jetTruthFinder(jetconstituents, ReconstructedParticles, Particle)")
    df = df.Define("jets_truthv2", "FCCAnalyses::jetTruthFinderV2(jet_p4, Particle)")

    #MC data

    #Truth Higgs collection
    df = df.Define("truth_H","FCCAnalyses::MCParticle::sel_pdgID(25,true)(Particle)")

    #Truth C quark collection 
    df = df.Define("truth_C","FCCAnalyses::MCParticle::sel_pdgID(4,true)(Particle)")

    #Truth B quark collection
    df = df.Define("truth_B","FCCAnalyses::MCParticle::sel_pdgID(5,true)(Particle)")
    

    p_vars = ["e","p","pt","px","py","pz","phi","mass","eta","tlv"]
    
    for p_var in p_vars:
        df = df.Define("truth_H_{}".format(p_var), "FCCAnalyses::MCParticle::get_{}(truth_H)".format(p_vars))

        df = df.Define("truth_C_{}".format(p_var), "FCCAnalyses::MCParticle::get_{}(truth_C)".format(p_vars))
        
        df = df.Define("truth_B_{}".format(p_var), "FCCAnalyses::MCParticle::get_{}(truth_B)".format(p_vars))


    #Retrieve H energy from tlv
    df = df.Define("truth_H_e_tlv", "FCCAnalyses::TLVHelpers::get_e(truth_H_tlv)")


    #MC Charm quark Data

    #MC Bottom quark Data

    return df


# Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis:
    # __________________________________________________________
    # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df = analysis_sequence(df)
        df = jet_sequence(df, rad, alg, sort, ecut)
        #df = jet_sequence(df, njets, exclusive) # again, was playing with exclusive parameter here. 
        # Don't remember if you need to pass it here.
        

        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = []

        # jets
        branches_jet = list(variables_jet.keys())
        branchList = branches_jet 

        branchList += antiktFlavourHelper.outputBranches()

        branchList += ["event_njet"]
        
        # branchList += ["all_invariant_masses_ee_kt"]
        branchList += ["all_invariant_masses"]

        branchList += ["recojetpair_isC"]
        branchList += ["recojetpair_isB"]

        # branchList += ["recoil_masses_ee_kt"]
        branchList += ["recoil_masses"]

        branchList += ["jet_e_corr"]
        branchList += ["jet_px_corr"]
        branchList += ["jet_py_corr"]
        branchList += ["jet_pz_corr"]
        
        
        # not corrected pt, e
        branchList += ["recojet_e"]
        branchList += ["recojet_px"]
        branchList += ["recojet_py"]
        branchList += ["recojet_pz"]
        
    
         # truth info
        branchList += ["jets_truth"]
        branchList += ["jets_truthv2"]
        
        # truth vars 
        truth_vars = ["e","p","pt","px","py","pz","phi","mass","eta"]
        for truth_var in truth_vars:
            branchList += [f"truth_H_{truth_var}"]
            branchList += [f"truth_C_{truth_var}"]
            branchList += [f"truth_B_{truth_var}"]

        branchList += "truth_H_e_tlv"

        #branchList += ["truth_H_e"]
        # branchList += ["truth_H_pt"]
        # branchList += ["truth_H_px"]
        # branchList += ["truth_H_py"]
        # branchList += ["truth_H_pz"]
        # branchList += ["truth_H_phi"]
        # branchList += ["truth_H_mass"]
        # branchList += ["truth_H_eta"]
        # branchList += ["truth_H_tlv"]  


        # vis kinematics
        branchList += ["vis_theta"]
        branchList += ["vis_M"]
        branchList += ["vis_E"]

      
        # leptons
        branchList += ["event_nel"]
        branchList += ["event_nmu"]
        branchList += ["muons_p"]
        branchList += ["electrons_p"]

        # MET
        MET_vars = ["e", "p", "pt", "px", "pt", "pz", "eta", "theta", "phi"]
        for MET_var in MET_vars:
            branchList += [f"RecoMissingEnergy_{MET_var}"]

        branchList = sorted(list(set(branchList))) # remove duplicates, sort

        return branchList
