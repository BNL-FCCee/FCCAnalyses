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
    JobName = values["JobName"]

print("batch:",batch)
print("EOSoutput:",EOSoutput)
print("name:",JobName)

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

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/" 

if(EOSoutput):
    outputDir = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/stage1/"
    eosType = "eosuser"
else:
    outputDir   = f"{JobName}/stage1/"

nCPUS       = 4
runBatch    = batch
batchQueue = "longlunch"
#compGroup = "group_u_FCC.local_gen"




""" 
Function to return all di-jet candidates, not just one near a certain resonance
"""

# resonanceBuilder::resonanceBuilder(float arg_resonance_mass, bool arg_return_all) {m_resonance_mass = arg_resonance_mass; return_all = arg_return_all;}
# ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> resonanceBuilder::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs) {


import ROOT
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> dijet_pair_all(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs){
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;

  int n = legs.size();
  if (n >1) {
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      edm4hep::ReconstructedParticleData reso;
      TLorentzVector reso_lv;
      for (int i = 0; i < n; ++i) {
          if (v[i]) {
            reso.charge += legs[i].charge;
            TLorentzVector leg_lv;
            leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
            reso_lv += leg_lv;
          }
      }
      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);
    } while (std::next_permutation(v.begin(), v.end()));
  }

  return result;

}

""") ###
#

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
        df2 = (
            df
            # define an alias for jet index collection
            # only Jet#2 and Jet#3 have something with non-negative collectionIDs...
            #.Alias("JetCollection", "Jet#2.index") # Average ~ 40 per event for ZccHbb
            #.Alias("JetCollection", "Jet#3.index")  # Average 0-4 per event for ZccHmumu...
            # "Jet" branch size is about 2 or 3 on average for ZccHmumu. Makes sense. Matches size of Jet#3 branch.
            # "Jet#2" branch size is about 40-50 on average....does not seem there

            # https://github.com/HEP-FCC/FCCAnalyses/blob/0aa81c44d7811d30fc6b098b8519b919461326d8/examples/FCCee/top/hadronic/analysis_stage1.py
            # https://github.com/HEP-FCC/FCCAnalyses/blob/7b02bc788ef14eb21201ac3bff2b67870208b29d/examples/basics/read_EDM4HEP.py

            # all particles 
            #.Alias("AllRecoParticles", "")

            #.Alias("GenParticles", "")

            # define the jet collection 
            #.Define("jets", "ReconstructedParticle::get(JetCollection, ReconstructedParticles)")

            .Define("n_jets", "ReconstructedParticle::get_n(Jet)") #count how many jets are in the event in total

            # Apparently "Jet" is an accesible collection.

            # define collections for jet flavour, and mother particle? 

            # select jets on pT 
            #.Define("selected_jets", "ReconstructedParticle::sel_pt(1)(jets)") # Stay loose with selection at this stage
            .Define("selected_jets", "ReconstructedParticle::sel_pt(1)(Jet)") # Loosest selection at this stage

            # create branches with jet (transverse momentum, rapidity, total momentum, energy)
            .Define("selected_jets_pt", "ReconstructedParticle::get_pt(selected_jets)")
            .Define("selected_jets_y", "ReconstructedParticle::get_y(selected_jets)")
            .Define("selected_jets_p", "ReconstructedParticle::get_p(selected_jets)")
            .Define("selected_jets_e", "ReconstructedParticle::get_e(selected_jets)")

            # create branches with jet (transverse momentum, rapidity, total momentum, energy)
            .Define("selected_jets_pt_0", "selected_jets_pt[0]") # what does [0] take? leading? 
            .Define("selected_jets_y_0", "selected_jets_y[0]")
            .Define("selected_jets_p_0", "selected_jets_p[0]")
            .Define("selected_jets_e_0", "selected_jets_e[0]")

            # find zed candidates from di-jet resonances

            .Define("dijet_pair_all", "dijet_pair_all(selected_jets)")

            #.Define("dijet_pair", "ReconstructedParticle::resonanceBuilder(91, 1)(selected_jets)") # seems to return one dijet pair. 
            #.Define("dijet_pair_nearZpeak", "ReconstructedParticle::resonanceBuilder(91, 0)(selected_jets)") # seems to return one dijet pair. 
            .Define("dijet_pair_nearZpeak", "ReconstructedParticle::resonanceBuilder(91)(selected_jets)") # seems to return one dijet pair. 

            # create branches with dijet (mass, transverse momentum, recoil, recoil mass, charge)
            .Define("dijet_pair_all_m", "ReconstructedParticle::get_mass(dijet_pair_all)")
            .Define("dijet_pair_all_pt", "ReconstructedParticle::get_pt(dijet_pair_all)")
            .Define("dijet_pair_all_recoil", "ReconstructedParticle::recoilBuilder(240)(dijet_pair_all)")
            .Define("dijet_pair_all_recoil_m", "ReconstructedParticle::get_mass(dijet_pair_all_recoil)")

            # create branches with dijet (mass, transverse momentum, recoil, recoil mass, charge)
            .Define("dijet_pair_nearZpeak_m", "ReconstructedParticle::get_mass(dijet_pair_nearZpeak)")
            .Define("dijet_pair_nearZpeak_pt", "ReconstructedParticle::get_pt(dijet_pair_nearZpeak)")
            .Define("dijet_pair_nearZpeak_recoil", "ReconstructedParticle::recoilBuilder(240)(dijet_pair_nearZpeak)")
            .Define("dijet_pair_nearZpeak_recoil_m", "ReconstructedParticle::get_mass(dijet_pair_nearZpeak_recoil)")            

            # save number of dijet pairs, number of jets.
            .Define("N_selected_jets", "selected_jets_pt.size()")
            .Define("N_dijet_pair", "dijet_pair_all_m.size()") # size of mass vector

            # Filter at least one candidate 
            #.Filter("dijet_pair_recoil_m.size()>0")


            #Define a variable based on a custom filter
            #.Define("MyFilter", "myFilter(zed_leptonic_m)") # maybe choose flavour here.             

        )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [

            # jets 
            "N_selected_jets",
            "selected_jets_pt",
            "selected_jets_y",
            "selected_jets_p",
            "selected_jets_e",

            # dijets pairs 
            "N_dijet_pair",
            "dijet_pair_all_pt",
            "dijet_pair_all_m",
            "dijet_pair_all_recoil_m",

            # Zed candidate near z mass peak 
            "dijet_pair_nearZpeak_m",           
            "dijet_pair_nearZpeak_pt",           
            "dijet_pair_nearZpeak_recoil_m",     
        ]
        return branchList
