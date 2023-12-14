"""
January 2023
Abraham Tishelman-Charny

The purpose of this python module is to perform initial selections and variable definitions for processing FCC files.
"""

import os
import urllib.request
import yaml

configFile = "RunConfig.yaml" # for the moment, need to specify full path so that HTCondor node can find this file (since afs is mounted). Need to check how to pass this as an input file to HTCondor job.
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
    #'wzp6_ee_ccH_HWW_ecm240':{'chunks':20},
    #'wzp6_ee_ccH_Hgg_ecm240' : {'chunks':20},
    #'wzp6_ee_ccH_HZa_ecm240' : {'chunks':20},
    #'wzp6_ee_ccH_Hss_ecm240' : {'chunks':20},
    #'wzp6_ee_ccH_Hmumu_ecm240':{'chunks':20},
    #'wzp6_ee_ccH_HZZ_ecm240' : {'chunks':20},	
    #'wzp6_ee_ccH_Htautau_ecm240' : {'chunks':20},
    #'wzp6_ee_ccH_Haa_ecm240' : {'chunks':20},
    #'wzp6_ee_ccH_Hcc_ecm240' : {'chunks':20},
    'wzp6_ee_ccH_Hbb_ecm240':{'chunks':20},

    # backgrounds. Option: 'fraction' : frac_value
    #'p8_ee_WW_ecm240' : {'chunks':3740},
    #'p8_ee_ZZ_ecm240' : {'chunks':562},
    #'p8_ee_Zqq_ecm240' : {'chunks':1007}
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/" 
procDict = "FCCee_procDict_winter2023_IDEA.json" 

if(EOSoutput):
    outputDir = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/stage1/"
    #outputDirEos = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/stage1/"
    #eosType = "eosuser"
else:
    outputDir   = f"{JobName}/stage1/"

nCPUS       = 4
runBatch    = batch
batchQueue = "workday" 

# Define any functionality which is not implemented in FCCAnalyses

print("test")

import ROOT
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<double> SumFlavorScores(ROOT::VecOps::RVec<double> recojet_isFlavor) {

    double score_1, score_2, pair_score; 
    ROOT::VecOps::RVec<double> recojetpair_isFlavor;

    // cannot compute any mass pair flavour score values, return a single non-physical value
    if(recojet_isFlavor.size() < 2){
        recojetpair_isFlavor.push_back(-99);
        return recojetpair_isFlavor; 
    }


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



ROOT::VecOps::RVec<double> all_recoil_masses(ROOT::VecOps::RVec<TLorentzVector> all_jet_4vectors){
  
  double m_sqrts = 240;
  auto recoil_p4 = TLorentzVector(0, 0, 0, m_sqrts);
  TLorentzVector tv1, tv2, tvpair; 
  double E, px, py, pz, recoil_mass;
  ROOT::VecOps::RVec<double> recoil_masses;

  // cannot compute any mass pair values, return a single non-physical value
  if(all_jet_4vectors.size() < 2){
    recoil_masses.push_back(-99);
    return recoil_masses;  
  }

    // For each jet, take its recoil mass using the remaining jets. Stop at last jet.
    for(int i = 0; i < all_jet_4vectors.size()-1; ++i) {

        tv1 = all_jet_4vectors.at(i);

        for(int j=i+1; j < all_jet_4vectors.size(); ++j){ // go until end

            tv2 = all_jet_4vectors.at(j); 
            E = tv1.E() + tv2.E();
            px = tv1.Px() + tv2.Px();
            py = tv1.Py() + tv2.Py();
            pz = tv1.Pz() + tv2.Pz();

            tvpair.SetPxPyPzE(px, py, pz, E);

            recoil_p4 = TLorentzVector(0, 0, 0, m_sqrts);
            recoil_p4 -= tvpair; 

            recoil_mass = recoil_p4.M();
            recoil_masses.push_back(recoil_mass);

        }
    }

  return recoil_masses;

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

## input file needed for unit test in CI
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

from FCCAnalyses.examples.FCCee.weaver.config import (
    variables_pfcand,
    variables_jet,
    variables_event, # assumes at least 2 jets for event_invariant_mass variable
)

from FCCAnalyses.addons.ONNXRuntime.python.jetFlavourHelper import JetFlavourHelper
from FCCAnalyses.addons.FastJet.python.jetClusteringHelper import ExclusiveJetClusteringHelper

jetFlavourHelper = None
jetClusteringHelper = None

print("tick 2")

# Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis:
    # __________________________________________________________
    # Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        global jetClusteringHelper
        global jetFlavourHelper

        from FCCAnalyses.examples.FCCee.weaver.config import collections

        tag = ""

        ## define jet clustering parameters
        jetClusteringHelper = ExclusiveJetClusteringHelper(collections["PFParticles"], njets, tag)

        ## run jet clustering
        df = jetClusteringHelper.define(df)

        print("tick 3")
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
        

#        df = df.Define("constituents", "JetConstituentsUtils::build_constituents_cluster({}, {})".format(jetClusteringHelper.constituents, 1))
#        df = df.Define("constituents_p4", "JetConstituentsUtils::get_constituents({}, {})".format(jetClusteringHelper.constituents, [0,1,2,3]))
        df = df.Define("all_invariant_masses", "JetConstituentsUtils::all_invariant_masses(jet_p4)")
        df = df.Define("recoil_masses", "all_recoil_masses(jet_p4)")
        df = df.Define("input_coll", "{}".format(jetClusteringHelper.input_coll))

        df = df.Define("MC_PDG", "FCCAnalyses::MCParticle::get_pdg(Particle)")
        df = df.Define("MC_genStatus", "FCCAnalyses::MCParticle::get_genStatus(Particle)")


        ## tagger inference
        df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df) 

        ## define variables using tagger inference outputs
#        df = df.Define("recojetpair_isC", "SumFlavorScores(recojet_isC)") 
#        df = df.Define("recojetpair_isB", "SumFlavorScores(recojet_isB)") 


        return df

    # __________________________________________________________
    # Mandatory: output function, please make sure you return the branchlist as a python list
    def output():

        branchList = []

        branches_jet = list(variables_jet.keys())
        branches_event = list(variables_event.keys())
        branches_pfcand = list(variables_pfcand.keys()) # extra info 

        branchList = branches_event  + branches_jet
        branchList += jetFlavourHelper.outputBranches()
        branchList += jetClusteringHelper.outputBranches()
        branchList += branches_pfcand
        branchList += ["all_invariant_masses"]
#       	branchList += ["recojetpair_isC"]
#        branchList += ["recojetpair_isB"]
        branchList += ["recoil_masses"]
        branchList += ["jet_p4"]
        branchList += ["input_coll"]
        branchList += ["MC_PDG"]
        branchList += ["MC_genStatus"]
#        branchList += ["constituents_1"]

        # remove duplicates 
        branchList = list(set(branchList))

        return branchList
