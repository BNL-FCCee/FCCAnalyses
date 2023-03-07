import os
import urllib.request
import yaml 

configFile = "ZccH/RunConfig.yaml"
with open(configFile, 'r') as cfg:
    values = yaml.safe_load(cfg)
    
    batch = values["batch"]
    EOSoutput = values["EOSoutput"]
    JobName = values["JobName"]
    njets = values["njets"]

print("batch:",batch)
print("EOSoutput:",EOSoutput)
print("name:",JobName)
print("njets:",njets)

processList = {
    # Z(cc)H by higgs final state 
    # 'wzp6_ee_ccH_HWW_ecm240':{'chunks':20},
    # 'wzp6_ee_ccH_Hgg_ecm240' : {'chunks':20},
    # 'wzp6_ee_ccH_HZa_ecm240' : {'chunks':20},
    # 'wzp6_ee_ccH_Hss_ecm240' : {'chunks':20},
    #'wzp6_ee_ccH_Hcc_ecm240' : {'chunks':20},
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

if(EOSoutput):
    outputDir = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/stage1/"
    eosType = "eosuser"
else:
    outputDir   = f"{JobName}/stage1/"

nCPUS       = 64
runBatch    = batch
batchQueue = "longlunch"

# Define any functionality which is not implemented in FCCAnalyses

import ROOT
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<double> SumFlavorScores(ROOT::VecOps::RVec<double> recojet_isFlavor) {

    double score_1, score_2, pair_score; 
    ROOT::VecOps::RVec<double> recojetpair_isFlavor;

    // For each jet, take its flavor score sum with the remaining jets. Stop at last jet.
    for(int i = 0; i < recojet_isFlavor.size()-1; ++i) {

    score_1 = recojet_isFlavor.at(i); 

        for(int j=i+1; j < recojet_isFlavor.size(); ++j){ // go until end
            score_2 = recojet_isFlavor.at(j);
            pair_score = score_1 + score_2; 
            recojetpair_isFlavor.push_back(pair_score);

        }
    }

    return recojetpair_isFlavor;
}

""") 

# ____________________________________________________________
def get_file_path(url, filename):
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        urllib.request.urlretrieve(url, os.path.basename(url))
        return os.path.basename(url)

# ____________________________________________________________

## input file needed for unit test in CI
testFile = "https://fccsw.web.cern.ch/fccsw/testsamples/wzp6_ee_nunuH_Hss_ecm240.root"

## latest particle transformer model, trainied on 9M jets in winter2023 samples - need to separate train/test samples?
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

from examples.FCCee.weaver.config import (
    variables_pfcand,
    variables_jet,
    variables_event,
)

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

        from examples.FCCee.weaver.config import collections

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
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets({})".format(jetClusteringHelper.jets))
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])") 
        df = df.Define("all_invariant_masses", "JetConstituentsUtils::AllInvariantMasses(jet_p4)")

        ## tagger inference
        df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df) #.Define("recojetpair_isC", "SumFlavorScores(recojet_isC)")
        df = df.Define("recojetpair_isC", "SumFlavorScores(recojet_isC)") # want to compute sum of "isC" scores for each jet pair. Important that you compute this in the same order as jet invariant masses.
        df = df.Define("recojetpair_isB", "SumFlavorScores(recojet_isB)") # want to compute sum of "isB" scores for each jet pair. Important that you compute this in the same order as jet invariant masses.

        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():

        branches_pfcand = list(variables_pfcand.keys())
        branches_jet = list(variables_jet.keys())
        branches_event = list(variables_event.keys())

        branchList = branches_event + branches_jet + branches_pfcand
        branchList += jetFlavourHelper.outputBranches()

        branchList += ["recojetpair_isC"]
        branchList += ["recojetpair_isB"]

        return branchList
