#Mandatory: List of processes

processList = {# 'p8_ee_ZZ_ecm240':{}, #Run the full statistics in one output file named <outputDir>/p8_ee_ZZ_ecm240.root
              # 'p8_ee_WW_ecm240':{} 
            #'p8_ee_ZH_#ecm240':{'fraction':0.2, 'output':'p8_ee_ZH_ecm240_out'} #Run 20% of the statistics in one file named <outputDir>/p8_ee_ZH_ecm240_out.root (example on how to change the output name)
              'wzp6_ee_ccH_Hbb_ecm240':{}#:{'output':'wzp6_ee_ccH_Hbb_ecm240_out'},
            #'wzp6_ee_bbH_Hcc_ecm240':{'output':'wzp6_ee_bbH_Hcc_ecm240_out'}
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "/usatlas/atlas01/atlasdisk/users/ivelisce/winter2023/IDEA/" #IDEA concept de detecteur

#Optional: output directory, default is local running directory
#outputDir   = "outputs/fccee/higgs/mH-recoil/hzz/stage1/"
outputDir   = "./2DPlotting"


#Optional: analysisName, default is ""
analysisName = "Testanalysis"

#Optional: ncpus, default is 4
nCPUS= 128

#anti-kt Radius 
#rad = 0.4

#Optional test file
#testFile ="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/p8_ee_ZH_ecm240/events_101027117.root"

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df
            .Alias("Muon0", "Muon#0.index")
            .Alias("Electron0", "Electron#0.index")

            .Define("muons",                "ReconstructedParticle::get(Muon0, ReconstructedParticles)")
            .Define("electrons",            "ReconstructedParticle::get(Electron0, ReconstructedParticles)")

#choose isolated muons and electrons in a specific momentum interval - here between 25 and 80 GeV- these should not become jets
            .Define("selected_muons",       "ReconstructedParticle::sel_p(25,80)(muons)")
            .Define("selected_electrons",   "ReconstructedParticle::sel_p(25,80)(electrons)")
            .Define("selected_leptons",     "ReconstructedParticle::merge(selected_muons, selected_electrons)")

#pt of these particules
            
            .Define("selected_muons_pt",    "ReconstructedParticle::get_pt(selected_muons)")
            .Define("selected_electrons_pt","ReconstructedParticle::get_pt(selected_electrons)")
            .Define("selected_leptons_pt",  "ReconstructedParticle::get_pt(selected_leptons)")

#px 

            .Define("selected_muons_px",    "ReconstructedParticle::get_px(selected_muons)")
            .Define("selected_electrons_px","ReconstructedParticle::get_px(selected_electrons)")
            .Define("selected_leptons_px",  "ReconstructedParticle::get_px(selected_leptons)")

#py 

            .Define("selected_muons_py",    "ReconstructedParticle::get_py(selected_muons)")
            .Define("selected_electrons_py","ReconstructedParticle::get_py(selected_electrons)")
            .Define("selected_leptons_py",  "ReconstructedParticle::get_py(selected_leptons)")

#pz 

            .Define("selected_muons_pz",    "ReconstructedParticle::get_pz(selected_muons)")
            .Define("selected_electrons_pz","ReconstructedParticle::get_pz(selected_electrons)")
            .Define("selected_leptons_pz",  "ReconstructedParticle::get_pz(selected_leptons)")

#rapidity
            
            .Define("selected_muons_y",     "ReconstructedParticle::get_y(selected_muons)")
            .Define("selected_electrons_y", "ReconstructedParticle::get_y(selected_electrons)")
            .Define("selected_leptons_y",   "ReconstructedParticle::get_y(selected_leptons)")


#momentum
            
            .Define("selected_muons_p",     "ReconstructedParticle::get_p(selected_muons)")
            .Define("selected_electrons_p", "ReconstructedParticle::get_p(selected_electrons)")
            .Define("selected_leptons_p",   "ReconstructedParticle::get_p(selected_leptons)")


#energy

            .Define("selected_muons_e",     "ReconstructedParticle::get_e(selected_muons)")
            .Define("selected_electrons_e", "ReconstructedParticle::get_e(selected_electrons)")
            .Define("selected_leptons_e",   "ReconstructedParticle::get_e(selected_leptons)")


#number of muons and electrons and their sum (merge: leptons) that were selected

            .Define("N_selected_muons",     "ReconstructedParticle::get_n(selected_muons)")
            .Define("N_selected_electrons", "ReconstructedParticle::get_n(selected_electrons)")
            .Define("N_selected_leptons",   "ReconstructedParticle::get_n(selected_leptons)")


#determine if loose leptons are necessary to code
#select muons, electrons, leptons with momentum 10GeV (to maximum value?)
            .Define("LooseMuons",       "ReconstructedParticle::sel_p(10)(muons)")
            .Define("LooseElectrons",   "ReconstructedParticle::sel_p(10)(electrons)")
            .Define("LooseLeptons",     "ReconstructedParticle::merge(LooseMuons, LooseElectrons)")

#find muons, electrons, leptons with momentum 2GeV (to maximum value?)
            .Define("LooseMuons_2",       "ReconstructedParticle::sel_p(2)(muons)")                  
            .Define("LooseElectrons_2",   "ReconstructedParticle::sel_p(2)(electrons)")
            .Define("LooseLeptons_2",     "ReconstructedParticle::merge(LooseMuons_2, LooseElectrons_2)") 
	    
#find muons, electrons, leptons with momentum 1GeV (to maximum value?)
            .Define("LooseMuons_1",       "ReconstructedParticle::sel_p(1)(muons)")                  
            .Define("LooseElectrons_1",   "ReconstructedParticle::sel_p(1)(electrons)")
            .Define("LooseLeptons_1",     "ReconstructedParticle::merge(LooseMuons_1, LooseElectrons_1)") 

#find number of three kinds of "loose leptons"
	        .Define("N_LooseLeptons", "ReconstructedParticle::get_n(LooseLeptons)")
            .Define("N_LooseLeptons_2", "ReconstructedParticle::get_n(LooseLeptons_2)")
            .Define("N_LooseLeptons_1", "ReconstructedParticle::get_n(LooseLeptons_1)")

#find transverse momentum, theta, phi, and momentum of all leptons with momentum 10GeV (to maximum value?)
            .Define("LooseLeptons_pt", "ReconstructedParticle::get_pt(LooseLeptons)")
            .Define("LooseLeptons_theta", "ReconstructedParticle::get_theta(LooseLeptons)")
            .Define("LooseLeptons_phi", "ReconstructedParticle::get_phi(LooseLeptons)")
            .Define("LooseLeptons_p", "ReconstructedParticle::get_p(LooseLeptons)")


#we now select 2 muons (or 0 if there are no muons)that are the best candidates for the Z using findElectrons
            .Define("zed_muons",         "ReconstructedParticle::findZleptons(selected_muons)")
#we use resonsanceBuilder to create the Z/sum with the selected pair of muons
            .Define("zed_muonic",        "ReconstructedParticle::resonanceBuilder(91)(zed_muons)") 
#we then remove the pair of muons from the total selected muons so that we can reiterate the process until there are 3 pairs of muons from the 3 Z*
            .Define("sselected_muons",   "ReconstructedParticle::remove(selected_muons, zed_muons)")

#we repeat the process and find a new muon pair out of the remaining selected muons
            .Define("zed_muonsbis",      "ReconstructedParticle::findZleptons(sselected_muons)")
#create the z/sum
            .Define("zed_muonic2",       "ReconstructedParticle::resonanceBuilder(91)(zed_muonsbis)")
#then we remove the muon pair from the total selected muons
            .Define("ssselected_muons",  "ReconstructedParticle::remove(sselected_muons, zed_muonsbis)")

#we repeatthe process again 
            .Define("zed_muonsbisbis",   "ReconstructedParticle::findZleptons(ssselected_muons)")
            .Define("zed_muonic3",       "ReconstructedParticle::resonanceBuilder(91)(zed_muonsbisbis)")

#we make a list of the muons in the momentum range that did not form pairs 
#by removing muon pairs from the remaning selected muons
            .Define("extramuons",        "ReconstructedParticle::remove(ssselected_muons, zed_muonsbisbis)")
#find the number of extra muons
            .Define("N_extramuons",      "ReconstructedParticle::get_n(extramuons)")


#we do the same procedure but with the electrons
            .Define("zed_electrons",     "ReconstructedParticle::findZleptons(selected_electrons)")
            .Define("zed_electronic",    "ReconstructedParticle::resonanceBuilder(91)(zed_electrons)")
            .Define("sselected_electrons", "ReconstructedParticle::remove(selected_electrons, zed_electrons)")

            .Define("zed_electronsbis",      "ReconstructedParticle::findZleptons(sselected_electrons)")
            .Define("zed_electronic2",     "ReconstructedParticle::resonanceBuilder(91)(zed_electronsbis)")
            .Define("ssselected_electrons","ReconstructedParticle::remove(sselected_electrons, zed_electronsbis)")

            .Define("zed_electronsbisbis", "ReconstructedParticle::findZleptons(ssselected_electrons)")
            .Define("zed_electronic3",     "ReconstructedParticle::resonanceBuilder(91)(zed_electronsbisbis)")


#we make a list of the elecrrons in the momentum range that did not form pairs 
#by removing muon pairs from the remaning selected electrons
            .Define("extraelectrons",      "ReconstructedParticle::remove(ssselected_electrons, zed_electronsbisbis)")
#find the number of extra muons

            .Define("N_extraelectrons",    "ReconstructedParticle::get_n(extraelectrons)")

#we want regroup the 3 Z together so the 3 muonic first
#we merge the first two pairs of muons 
            .Define("mergemuonic1",         "ReconstructedParticle::merge(zed_muonic, zed_muonic2)")
#then we merge those two pairs with the third pair
            .Define("mergemuonic2",         "ReconstructedParticle::merge(mergemuonic1, zed_muonic3)")

#then the 3 electronic pairs (there will be a maximum 3/6 non empty since 3 Z max, or background)
#we merge the first two pairs of electrons
            .Define("mergeelectronic1",     "ReconstructedParticle::merge(zed_electronic, zed_electronic2)")
#then we those two pairs with the third pair
            .Define("mergeelectronic2",     "ReconstructedParticle::merge(mergeelectronic1, zed_electronic3)")

#and here we merge the 6 pairs of leptons to create zed_leptonic
            .Define("zed_leptonic",         "ReconstructedParticle::merge(mergemuonic2, mergeelectronic2)")

#we make the list of all the leptons that are selected in the momentum range but did not form pairs 
            .Define("extraleptons", "ReconstructedParticle::merge(extramuons, extraelectrons)")

#then we find the number of elements in zed_leptonic Z and extraleptons
            .Define("N_zed_leptonic",       "ReconstructedParticle::get_n(zed_leptonic)")
            .Define("N_extraleptons",       "ReconstructedParticle::get_n(extraleptons)")
            

#properties of Z pairs
           #get_e does not work on zed_leptonic so we can do it on the leptons directly and then sum [0] and [1] to have the 
           #energy component of the threshold limit values (tlv) of the Z (we only have one candidate for the Z in stage 2 with the filter
           #so we do this with the first and best pairs of electrons and muons of findZleptons)

            #energy of best electron pair
            .Define("zed_electrons_e", "ReconstructedParticle::get_e(zed_electrons)")
            #energy of best muon pair
            .Define("zed_muons_e", "ReconstructedParticle::get_e(zed_muons)")
            .Define("zed_leptons", "ReconstructedParticle::merge(zed_electrons, zed_muons)")
        #energy of the best muon and electron pairs
            .Define("zed_leptons_e", "ReconstructedParticle::get_e(zed_leptons)")
            
            #mass of all candidate lepton pairs
            .Define("zed_leptonic_m",       "ReconstructedParticle::get_mass(zed_leptonic)")
            #same for pt
            .Define("zed_leptonic_pt",      "ReconstructedParticle::get_pt(zed_leptonic)")
            #same for px
            .Define("zed_leptonic_px",      "ReconstructedParticle::get_px(zed_leptonic)")
            #same for py
            .Define("zed_leptonic_py",      "ReconstructedParticle::get_py(zed_leptonic)") 
            #same for pz
            .Define("zed_leptonic_pz",      "ReconstructedParticle::get_pz(zed_leptonic)")
            #same for p
            .Define("zed_leptonic_p",       "ReconstructedParticle::get_p(zed_leptonic)")
            #recoil for all candidate lepton pairs
            .Define("zed_leptonic_recoil",  "ReconstructedParticle::recoilBuilder(240)(zed_leptonic)")
            #mass of recoil for all candidate lepton pairs
            .Define("zed_leptonic_recoil_m","ReconstructedParticle::get_mass(zed_leptonic_recoil)")  
            #charge of all lepton pairs          
            .Define("zed_leptonic_charge",  "ReconstructedParticle::get_charge(zed_leptonic)")
            #same for theta
            .Define("zed_leptonic_theta",   "ReconstructedParticle::get_theta(zed_leptonic)")
            #same for phi
            .Define("zed_leptonic_phi",     "ReconstructedParticle::get_phi(zed_leptonic)")
            #same for rapidity (y)
            .Define("zed_leptonic_y",       "ReconstructedParticle::get_y(zed_leptonic)")
            #same for eta
            .Define("zed_leptonic_eta",     "ReconstructedParticle::get_eta(zed_leptonic)")
            #same for cos 
            .Define("zed_leptonic_cos",     "cos(ReconstructedParticle::get_theta(zed_leptonic))")

    #we make the list of leptons for the pairs: so the selected - the extras
            .Define("taken_leptons",    "ReconstructedParticle::remove(selected_leptons, extraleptons)")
            .Define("N_taken_leptons",  "ReconstructedParticle::get_n(taken_leptons)")

    #we select all the particles used for the pairs of Z to reconstruct the jets
            .Define("my_recoparticles",  "ReconstructedParticle::remove(ReconstructedParticles, taken_leptons)")
    #find px reco particles (RP) 
            .Define("RP_px", "ReconstructedParticle::get_px(my_recoparticles)")
    #same for py
            .Define("RP_py", "ReconstructedParticle::get_py(my_recoparticles)")
    #same for pz
            .Define("RP_pz", "ReconstructedParticle::get_pz(my_recoparticles)")
    #find energy for RP
            .Define("RP_e", "ReconstructedParticle::get_e(my_recoparticles)")

#we create pseudo jets with the reco particles using set_pseudoJets method with px, py, pz, and, energy of RP
            .Define("pseudo_jets",  "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")


#the truth flavor reconstructed by the jets - be careful of the definition of get_flavour 
#the truth flavor is the falvor of the particle that gave birth to the jet that we know via pdgid MONTE CARLO 
#when we reconstruct the jets we cannot know which exactly is the origine of the jet so we have to determine the best via get_flavor which is not
#tagging but it is a way to attribute a truth value of flavor to a jet and not a particle which is not evident

    #Durham Algorithm for n=4
    #enter number of jets in second spot in parenthesis

     #exclusive Durham clustering with number of jets sets to 4
           # .Define("FCCAnalysesJets_ee_genkt4",  "JetClustering::clustering_ee_kt(2, 4, 1, 0)(pseudo_jets)")
         .Define("FCCAnalysesJets_ee_antikt4", "JetClustering::clustering_antikt(0.4,0,0,0,0)(pseudo_jets)")
         .Define("jets_ee_genkt4",  "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt4)")
    #get properties of jets 
            .Define("N_jets_4", "JetClusteringUtils::get_n(jets_ee_genkt4)")
            .Define("jets_px4",  "JetClusteringUtils::get_px(jets_ee_genkt4)")
            .Define("jets_py4",  "JetClusteringUtils::get_py(jets_ee_genkt4)")
            .Define("jets_pz4",  "JetClusteringUtils::get_pz(jets_ee_genkt4)")
            .Define("jets_p4",   "JetClusteringUtils::get_p(jets_ee_genkt4)")
            .Define("jets_e4",   "JetClusteringUtils::get_e(jets_ee_genkt4)")
            .Define("jets_m4",   "JetClusteringUtils::get_m(jets_ee_genkt4)")
            .Define("jets_pt4",  "JetClusteringUtils::get_pt(jets_ee_genkt4)")
            .Define("jets_y4",   "JetClusteringUtils::get_y(jets_ee_genkt4)")
            .Define("jets_eta4", "JetClusteringUtils::get_eta(jets_ee_genkt4)")
            .Define("jets_theta4", "JetClusteringUtils::get_theta(jets_ee_genkt4)")
            .Define("jets_phi4", "JetClusteringUtils::get_phi(jets_ee_genkt4)")

#get constituents of jets
            .Define("jetconstituents_ee_genkt4", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt4)")
#build constituent clusters using reco particles and jet constituents
            .Define("jetconstituents_ee_4", "JetConstituentsUtils::build_constituents_cluster(my_recoparticles, jetconstituents_ee_genkt4)"
)
#count constituents 
            .Define("jetconstituents_4", "JetConstituentsUtils::count_consts(jetconstituents_ee_4)")
#theta
            .Define("jetconstituents_4_theta", "JetConstituentsUtils::get_theta(jetconstituents_ee_4)")
#phi
            .Define("jetconstituents_4_phi", "JetConstituentsUtils::get_phi(jetconstituents_ee_4)")
#energy
            .Define("jetconstituents_4_energy", "JetConstituentsUtils::get_e(jetconstituents_ee_4)")

#dmerge for each of the jets
            .Define("dmerge_4_45", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt4, 4)")
            .Define("dmerge_4_34", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt4, 3)")
            .Define("dmerge_4_23", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt4, 2)")
            .Define("dmerge_4_12", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt4, 1)")

#flavour of jet
           .Define("jets_ee_genkt_flavour4",   "JetTaggingUtils::get_flavour(jets_ee_genkt4, Particle)")   
            .Define("jets_ee_flavour4", "JetTaggingUtils::get_flavour(jets_ee_genkt4, Particle)")
        
            #anti-kt Algorithm (currently R=0.4)
            #enter radius value R
            .Define("FCCAnalysesJets_ee_antikt4", "JetClustering::clustering_antikt(0.4,0,0,0,0)(pseudo_jets)")
            .Define("jets_antikt4", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_antikt4)")
            .Define("N_jets_antikt4", "JetClusteringUtils::get_n(jets_antikt4)")
            .Define("jets_antikt_e4", "JetClusteringUtils::get_e(jets_antikt4)")
            .Define("jets_antikt_px4", "JetClusteringUtils::get_px(jets_antikt4)")
            .Define("jets_antikt_py4", "JetClusteringUtils::get_py(jets_antikt4)")
            .Define("jets_antikt_pz4", "JetClusteringUtils::get_pz(jets_antikt4)")
            .Define("jets_antikt_theta4", "JetClusteringUtils::get_theta(jets_antikt4)")
            .Define("jets_antikt_phi4", "JetClusteringUtils::get_phi(jets_antikt4)")

            .Define("jetconstituents_ee_antikt4", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_antikt4)")
            .Define("jetconstituents_ee_antikt4_utile", "JetConstituentsUtils::build_constituents_cluster(my_recoparticles, jetconstituents_ee_antikt4)")
            .Define("jetconstituents_antikt4", "JetConstituentsUtils::count_consts(jetconstituents_ee_antikt4_utile)")
            .Define("jetconstituents_antikt4_theta", "JetConstituentsUtils::get_theta(jetconstituents_ee_antikt4_utile)")
            .Define("jetconstituents_antikt4_phi", "JetConstituentsUtils::get_phi(jetconstituents_ee_antikt4_utile)")
            .Define("jetconstituents_antikt4_energy", "JetConstituentsUtils::get_e(jetconstituents_ee_antikt4_utile)")

            .Define("dmerge_antikt4_45", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4, 4)")
            .Define("dmerge_antikt4_34", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4, 3)")
            .Define("dmerge_antikt4_23", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4, 2)")
            .Define("dmerge_antikt4_12", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4, 1)")
            #end anti-kt 
           
#hzz monte carlo
        
            .Alias("Particle1", "Particle#1.index")

            .Define("ZH_decay", "MCParticle::fill_ZH_decay(Particle, Particle1)")
            #select truth Z
            .Define("truth_Z", "MCParticle::sel_pdgID(23, true)(Particle)")
            #select truth H
            .Define("truth_H", "MCParticle::sel_pdgID(25, true)(Particle)")
            #truth Z theta
            .Define("truth_Z_theta", "MCParticle::get_theta(truth_Z)")
            #truth Z phi
            .Define("truth_Z_phi", "MCParticle::get_phi(truth_Z)")
            #truth H theta
            .Define("truth_H_theta", "MCParticle::get_theta(truth_H)")
            #truth H phi
            .Define("truth_H_phi", "MCParticle::get_phi(truth_H)")
        
            )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            "selected_muons_pt",
            "selected_electrons_pt",
            "selected_leptons_pt",

            "selected_muons_px",
            "selected_electrons_px",
            "selected_leptons_px",

            "selected_muons_py",
            "selected_electrons_py",
            "selected_leptons_py",

            "selected_muons_pz",
            "selected_electrons_pz",
            "selected_leptons_pz",

            "selected_muons_y",
            "selected_electrons_y",
            "selected_leptons_y",

            "selected_muons_p",
            "selected_electrons_p",
            "selected_leptons_p",

            "selected_muons_e",
            "selected_electrons_e",
            "selected_leptons_e",

            "N_selected_muons",
            "N_selected_electrons",
            "N_selected_leptons",
            "N_extramuons",
            "N_extraelectrons",
            "N_taken_leptons",
            "N_zed_leptonic",

            "zed_leptonic_pt",
            "zed_leptonic_px",
            "zed_leptonic_py", 
            "zed_leptonic_pz",
            "zed_leptonic_p",
            "zed_leptonic_m",
            "zed_leptonic_charge",
            "zed_leptonic_recoil_m",
            "zed_leptonic_phi",
            "zed_leptonic_theta",
            "zed_leptonic_cos",
            "zed_leptonic_y",
            "zed_leptonic_eta",
           
            "jets_px4",
            "jets_py4",
            "jets_pz4",
            "jets_p4",
            "jets_e4",
            "jets_m4",
            "jets_pt4",
            "jets_y4",
            "jets_eta4",
            "jets_theta4",
            "jets_phi4",
            "jetconstituents_ee_genkt4",
            "jetconstituents_ee_4",
            "jetconstituents_4",
            "dmerge_4_45",
            "dmerge_4_34",
            "dmerge_4_23",
            "dmerge_4_12",
            "jetconstituents_4_theta",
            "jetconstituents_4_phi",
            "jetconstituents_4_energy",
            "jets_ee_genkt_flavour4",
            "jets_ee_flavour4",

            "ZH_decay",
            
            "truth_Z_theta",
            "truth_Z_phi",
            
            "truth_H_theta",
            "truth_H_phi",
        
            "N_jets_antikt4",
            "jets_antikt_e4",
            "jets_antikt_px4",
            "jets_antikt_py4",
            "jets_antikt_pz4",
            "jets_antikt_theta4",
            "jets_antikt_phi4",
            "jetconstituents_antikt4",
            "jetconstituents_antikt4_theta",
            "jetconstituents_antikt4_phi",
            "jetconstituents_antikt4_energy",
            "dmerge_antikt4_45",
            "dmerge_antikt4_34",
            "dmerge_antikt4_23",
            "dmerge_antikt4_12",

            "N_jets_4",

            "N_LooseLeptons",
	        "N_LooseLeptons_2",
            "N_LooseLeptons_1",
            "LooseLeptons_pt",
            "LooseLeptons_theta", 
            "LooseLeptons_phi", 
            "LooseLeptons_p",
            
            "zed_electrons_e",
            "zed_muons_e",
            "zed_leptons_e"
            
]

        return branchList