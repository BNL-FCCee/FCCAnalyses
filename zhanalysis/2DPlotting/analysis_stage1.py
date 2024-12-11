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

#Number of jets (Durham kt)

#Radius for anti-kt
#R = 0.4

#Clustering Types: 
#algs = {"antikt", "cambridge","eekt", "ee_genkt"}

#choose which algorithm to use
#1=anti-kt, 2=Cambridge, 3=Durham kt, 4= kt
#alg=2


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

#we make the list of all the leptons that are selected in the momentum range but did not form pairs 
            .Define("extraleptons", "ReconstructedParticle::merge(extramuons, extraelectrons)")

    #we make the list of leptons for the pairs: so the selected - the extras
            .Define("taken_leptons",    "ReconstructedParticle::remove(selected_leptons, extraleptons)")

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



        #collection of reconstructed particles
            .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
        #collection of Monte Carlo particles
            .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")
            .Define('RP_MC_index',"ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles)")
            .Define()
            


    #enter number of jets in second spot in parenthesis

#exclusive Durham clustering with n=4
            .Define("FCCAnalysesJets_kt4",  "JetClustering::clustering_ee_kt(2, 4, 1, 0)(pseudo_jets)")
            .Define("jets_kt4",  "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_kt4)")

         #get properties of jets 
            .Define("jets_kt_px4",  "JetClusteringUtils::get_px(jets_kt4)")
            .Define("jets_kt_py4",  "JetClusteringUtils::get_py(jets_kt4)")
            .Define("jets_kt_pz4",  "JetClusteringUtils::get_pz(jets_kt4)")
            .Define("jets_kt_p4",   "JetClusteringUtils::get_p(jets_kt4)")
            .Define("jets_kt_e4",   "JetClusteringUtils::get_e(jets_kt4)")
            .Define("jets_kt_m4",   "JetClusteringUtils::get_m(jets_kt4)")
            .Define("jets_kt_pt4",  "JetClusteringUtils::get_pt(jets_kt4)")
            .Define("jets_kt_y4",   "JetClusteringUtils::get_y(jets_kt4)")
            .Define("jets_kt_eta4", "JetClusteringUtils::get_eta(jets_kt4)")
            .Define("jets_kt_theta4", "JetClusteringUtils::get_theta(jets_kt4)")
            .Define("jets_kt_phi4", "JetClusteringUtils::get_phi(jets_kt4)")

        #get number of constituents in each jet
            .Define("jetconstituents_get_kt4", "JetClusteringUtils::get_constituents(FCCAnalysesJets_kt4)")
        #build constituent clusters using reconstructed particles and number of jet constituents
            .Define("jetconstituents_kt_4", "JetConstituentsUtils::build_constituents_cluster(my_recoparticles, jetconstituents_get_kt4)")
            .Define("jetconstituents_kt4", "JetConstituentsUtils::count_consts(jetconstituents_kt_4)")

         #theta
            .Define("jetconstituents_kt4_theta", "JetConstituentsUtils::get_theta(jetconstituents_kt_4)")
         #phi
            .Define("jetconstituents_kt4_phi", "JetConstituentsUtils::get_phi(jetconstituents_kt_4)")
        #energy
            .Define("jetconstituents_kt4_energy", "JetConstituentsUtils::get_e(jetconstituents_kt_4)")


            #.Define("jets_alg_flavour4",   "JetTaggingUtils::get_flavour(jets_4, Particle)")   
            #.Define("jets_flavour4", "JetTaggingUtils::get_flavour(jets_4, Particle)")
            .Define("jets_kt_flavour_4",   "JetTaggingUtils::get_flavour(jets_kt4, Particle)")   
            .Define("jets_kt_flavour4", "JetTaggingUtils::get_flavour(jets_kt4, Particle)")
        

#antikt clustering with R=0.4
           .Define("FCCAnalysesJets_antikt4", "JetClustering::clustering_antikt(0.4,0,0,0,0)(pseudo_jets)")
           .Define("jets_antikt4",  "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_antikt4)")

        #get properties of jets 
            .Define("jets_antikt_px4",  "JetClusteringUtils::get_px(jets_antikt4)")
            .Define("jets_antikt_py4",  "JetClusteringUtils::get_py(jets_antikt4)")
            .Define("jets_antikt_pz4",  "JetClusteringUtils::get_pz(jets_antikt4)")
            .Define("jets_antikt_p4",   "JetClusteringUtils::get_p(jets_antikt4)")
            .Define("jets_antikt_e4",   "JetClusteringUtils::get_e(jets_antikt4)")
            .Define("jets_antikt_m4",   "JetClusteringUtils::get_m(jets_antikt4)")
            .Define("jets_antikt_pt4",  "JetClusteringUtils::get_pt(jets_antikt4)")
            .Define("jets_antikt_y4",   "JetClusteringUtils::get_y(jets_antikt4)")
            .Define("jets_antikt_eta4", "JetClusteringUtils::get_eta(jets_antikt4)")
            .Define("jets_antikt_theta4", "JetClusteringUtils::get_theta(jets_antikt4)")
            .Define("jets_antikt_phi4", "JetClusteringUtils::get_phi(jets_antikt4)")
        
        #get number of constituents in each jet
            .Define("jetconstituents_get_antikt4", "JetClusteringUtils::get_constituents(FCCAnalysesJets_antikt4)")
        #build constituent clusters using reconstructed particles and number of jet constituents
            .Define("jetconstituents_antikt_4", "JetConstituentsUtils::build_constituents_cluster(my_recoparticles, jetconstituents_get_antikt4)")

            .Define("jetconstituents_antikt4", "JetConstituentsUtils::count_consts(jetconstituents_antikt_4)")


        #theta
            .Define("jetconstituents_antikt4_theta", "JetConstituentsUtils::get_theta(jetconstituents_antikt_4)")
         #phi
            .Define("jetconstituents_antikt4_phi", "JetConstituentsUtils::get_phi(jetconstituents_antikt_4)")
        #energy
            .Define("jetconstituents_antikt4_energy", "JetConstituentsUtils::get_e(jetconstituents_antikt_4)")



        #the truth flavor reconstructed by the jets - be careful of the definition of get_flavour 
        #the truth flavor is the falvor of the particle that gave birth to the jet that we know via pdgid MONTE CARLO 

        #when we reconstruct the jets we cannot know which exactly is the origine of the jet so we have to determine the best via get_flavor which is not
        #tagging but it is a way to attribute a truth value of flavor to a jet and not a particle which is not evident
        #flavour of jet
            .Define("jets_antikt_flavour_4",   "JetTaggingUtils::get_flavour(jets_antikt4, Particle)")   
            .Define("jets_antikt_flavour4", "JetTaggingUtils::get_flavour(jets_antikt4, Particle)")
        
        

#hzz monte carlo
        
            .Alias("Particle1", "Particle#1.index")
            .Define("ZH_decay", "MCParticle::fill_ZH_decay(Particle, Particle1)")
            #select truth H
            .Define("truth_H", "MCParticle::sel_pdgID(25, true)(Particle)")
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
        #Durham Kt
            "jets_kt_px4",
            "jets_kt_py4",
            "jets_kt_pz4",
            "jets_kt_p4",
            "jets_kt_e4",
            "jets_kt_m4",
            "jets_kt_pt4",
            "jets_kt_y4",
            "jets_kt_eta4",
            "jets_kt_theta4",
            "jets_kt_phi4",

            "jetconstituents_get_kt4",
            "jetconstituents_kt_4",
            "jetconstituents_kt4",
            
            "jetconstituents_kt4_theta",
            "jetconstituents_kt4_phi",
            "jetconstituents_kt4_energy",
            
            "jets_kt_flavour_4",
            "jets_kt_flavour4",

        #Anti-kt R=0.4

            "jets_antikt_px4",
            "jets_antikt_py4",
            "jets_antikt_pz4",
            "jets_antikt_p4",
            "jets_antikt_e4",
            "jets_antikt_m4",
            "jets_antikt_pt4",
            "jets_antikt_y4",
            "jets_antikt_eta4",
            "jets_antikt_theta4",
            "jets_antikt_phi4",
        
            "jetconstituents_get_antikt4",
            "jetconstituents_antikt_4",
            "jetconstituents_antikt4",

            "jetconstituents_antikt4_theta",
            "jetconstituents_antikt4_phi",
            "jetconstituents_antikt4_energy",

            "jets_antikt_flavour_4",
            "jets_antikt_flavour4",
        ###
            "ZH_decay",
            "truth_H_theta",
            "truth_H_phi",

            
]

        return branchList