#Mandatory: List of processes

processList = {#'wzp6_ee_mumuH_HZZ_ecm240':{},
               #'wzp6_ee_mumuH_HWW_ecm240':{},
               #'wzp6_ee_mumuH_HZa_ecm240':{},
               #'wzp6_ee_mumuH_Haa_ecm240':{},
               #'wzp6_ee_mumuH_Hbb_ecm240':{},
               #'wzp6_ee_mumuH_Hcc_ecm240':{},
               #'wzp6_ee_mumuH_Hgg_ecm240':{},
               #'wzp6_ee_mumuH_Hmumu_ecm240':{},
               #'wzp6_ee_mumuH_Hss_ecm240':{},
               #'wzp6_ee_mumuH_Htautau_ecm240':{},
               #'wzp6_ee_mumu_ecm240':{},
               #'wzp6_ee_eeH_HZZ_ecm240':{},
               #'wzp6_ee_eeH_HWW_ecm240':{},
               #'wzp6_ee_eeH_HZa_ecm240':{},
               #'wzp6_ee_eeH_Haa_ecm240':{},
               #'wzp6_ee_eeH_Hbb_ecm240':{},
               #'wzp6_ee_eeH_Hcc_ecm240':{},
               #'wzp6_ee_eeH_Hgg_ecm240':{},
               #'wzp6_ee_eeH_Hmumu_ecm240':{},
               #'wzp6_ee_eeH_Hss_ecm240':{},
               #'wzp6_ee_eeH_Htautau_ecm240':{},
               #'wzp6_ee_ee_Mee_30_150_ecm240':{},
               'p8_ee_ZZ_ecm240':{}, #Run the full statistics in one output file named <outputDir>/p8_ee_ZZ_ecm240.root
               'p8_ee_WW_ecm240':{} #Run 50% of the statistics in two files named <outputDir>/p8_ee_WW_ecm240/chunk<N>.root
    #'p8_ee_ZH_#ecm240':{'fraction':0.2, 'output':'p8_ee_ZH_ecm240_out'} #Run 20% of the statistics in one file named <outputDir>/p8_ee_ZH_ecm240_out.root (example on how to change the output name)
               }

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/" #IDEA concept de detecteur

#Optional: output directory, default is local running directory
#outputDir   = "outputs/fccee/higgs/mH-recoil/hzz/stage1/"
# outputDir   = "outputs/fccee/higgs/mH-recoil/hzz/stage1/"
outputDir   = ""


#Optional: analysisName, default is ""
analysisName = "Testanalysis"

#Optional: ncpus, default is 4
nCPUS       = 128

#Optional running on HTCondor, default is False
#runBatch    = False

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"

#Optional test file
#testFile ="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/p8_ee_ZH_ecm240/events_101027117.root"

nominal_beff_WP80_b = 0.80
nominal_beff_WP80_c = 0.004
nominal_beff_WP80_l = 0.0005
nominal_beff_WP80_g = 0.007

#diapo13 lg-jettagging-fccee-krakow2023.pdf

nominal_ceff_WP80_b = 0.020
nominal_ceff_WP80_c = 0.80
nominal_ceff_WP80_l = 0.009
nominal_ceff_WP80_g = 0.025

#diapo13 lg-jettagging-fccee-krakow2023.pdf

nominal_geff_WP80_b = 0.020
nominal_geff_WP80_c = 0.050
nominal_geff_WP80_l = 0.150
nominal_geff_WP80_g = 0.80

#diapo14 de lg-jettagging-fccee-krakow2023.pdf

nominal_beff_WP90_b = 0.90
nominal_beff_WP90_c = 0.02
nominal_beff_WP90_l = 0.001
nominal_beff_WP90_g = 0.02

#diapo13 lg-jettagging-fccee-krakow2023.pdf

nominal_ceff_WP90_b = 0.040
nominal_ceff_WP90_c = 0.90
nominal_ceff_WP90_l = 0.070
nominal_ceff_WP90_g = 0.070

#diapo13 lg-jettagging-fccee-krakow2023.pdf

nominal_geff_WP90_b = 0.025
nominal_geff_WP90_c = 0.07
nominal_geff_WP90_l = 0.25
nominal_geff_WP90_g = 0.9

#diapo14 lg-jettagging-fccee-krakow2023.pdf

nominal_ceff_WP70_b = 0.009
nominal_ceff_WP70_c = 0.70
nominal_ceff_WP70_l = 0.002
nominal_ceff_WP70_g = 0.010


beff_b = nominal_beff_WP80_b
beff_c = nominal_beff_WP80_c
beff_l = nominal_beff_WP80_l
beff_g = nominal_beff_WP80_g

ceff_b = nominal_ceff_WP80_b
ceff_c = nominal_ceff_WP80_c
ceff_l = nominal_ceff_WP80_l
ceff_g = nominal_ceff_WP80_g

geff_b = nominal_geff_WP80_b
geff_c = nominal_geff_WP80_c
geff_l = nominal_geff_WP80_l
geff_g = nominal_geff_WP80_g

leff_l = 0.6
leff_b = 0.04
leff_c = 0.05
leff_g = 0.2

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df
            .Alias("Muon0", "Muon#0.index")
            .Alias("Electron0", "Electron#0.index")


#             .Define("vrai_Z", "MCParticle::sel_pdgID(23, true)(Particle)")
            .Define("vrai_Z", "MCParticle::sel_pdgID(23, true)(Particle)")
            .Define("vrai_H", "MCParticle::sel_pdgID(25, true)(Particle)")


            .Define("muons",                "ReconstructedParticle::get(Muon0, ReconstructedParticles)")
            .Define("electrons",            "ReconstructedParticle::get(Electron0, ReconstructedParticles)")

#on sélectionne des muons/électrons qui sont isolés + dans un certain intervalle de momentum : ne devraient pas venir des jets 

            .Define("selected_muons",       "ReconstructedParticle::sel_p(25,80)(muons)")
            .Define("selected_electrons",   "ReconstructedParticle::sel_p(25,80)(electrons)")
            .Define("selected_leptons",     "ReconstructedParticle::merge(selected_muons, selected_electrons)")

#pt de ces particules
            
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

#rapidité
            
            .Define("selected_muons_y",     "ReconstructedParticle::get_y(selected_muons)")
            .Define("selected_electrons_y", "ReconstructedParticle::get_y(selected_electrons)")
            .Define("selected_leptons_y",   "ReconstructedParticle::get_y(selected_leptons)")


#momentum
            
            .Define("selected_muons_p",     "ReconstructedParticle::get_p(selected_muons)")
            .Define("selected_electrons_p", "ReconstructedParticle::get_p(selected_electrons)")
            .Define("selected_leptons_p",   "ReconstructedParticle::get_p(selected_leptons)")


#energie

            .Define("selected_muons_e",     "ReconstructedParticle::get_e(selected_muons)")
            .Define("selected_electrons_e", "ReconstructedParticle::get_e(selected_electrons)")
            .Define("selected_leptons_e",   "ReconstructedParticle::get_e(selected_leptons)")


#nombre de muons/électrons et de la somme (merge : leptons) qui ont passé la sélection

            .Define("N_selected_muons",     "ReconstructedParticle::get_n(selected_muons)")
            .Define("N_selected_electrons", "ReconstructedParticle::get_n(selected_electrons)")
            .Define("N_selected_leptons",   "ReconstructedParticle::get_n(selected_leptons)")


            .Define("LooseMuons",       "ReconstructedParticle::sel_p(10)(muons)")
            .Define("LooseElectrons",   "ReconstructedParticle::sel_p(10)(electrons)")
            .Define("LooseLeptons",     "ReconstructedParticle::merge(LooseMuons, LooseElectrons)")

            .Define("LooseMuons_2",       "ReconstructedParticle::sel_p(2)(muons)")                  
            .Define("LooseElectrons_2",   "ReconstructedParticle::sel_p(2)(electrons)")
            .Define("LooseLeptons_2",     "ReconstructedParticle::merge(LooseMuons_2, LooseElectrons_2)") 
	    
            .Define("LooseMuons_1",       "ReconstructedParticle::sel_p(1)(muons)")                  
            .Define("LooseElectrons_1",   "ReconstructedParticle::sel_p(1)(electrons)")
            .Define("LooseLeptons_1",     "ReconstructedParticle::merge(LooseMuons_1, LooseElectrons_1)") 


	    .Define("N_LooseLeptons", "ReconstructedParticle::get_n(LooseLeptons)")
            .Define("N_LooseLeptons_2", "ReconstructedParticle::get_n(LooseLeptons_2)")
            .Define("N_LooseLeptons_1", "ReconstructedParticle::get_n(LooseLeptons_1)")

            .Define("LooseLeptons_pt", "ReconstructedParticle::get_pt(LooseLeptons)")
            .Define("LooseLeptons_theta", "ReconstructedParticle::get_theta(LooseLeptons)")
            .Define("LooseLeptons_phi", "ReconstructedParticle::get_phi(LooseLeptons)")
            .Define("LooseLeptons_p", "ReconstructedParticle::get_p(LooseLeptons)")

#on sélectionne 2 muons (ou 0 si pas de muons) étant les meilleurs candidats au Z ; findZleptons garde les leptons, ne construit pas le Z
#on crée le Z/somme avec resonanceBuilder
#on enlève ces deux muons des selected muons et on réitère la procédure puisque l'on peut avoir jusqu'à trois paires de muons venant des 3 Z*

            .Define("zed_muons",         "ReconstructedParticle::findZleptons(selected_muons)")
            .Define("zed_muonic",        "ReconstructedParticle::resonanceBuilder(91)(zed_muons)") 
            .Define("sselected_muons",   "ReconstructedParticle::remove(selected_muons, zed_muons)")
            .Define("zed_muonsbis",      "ReconstructedParticle::findZleptons(sselected_muons)")
            .Define("zed_muonic2",       "ReconstructedParticle::resonanceBuilder(91)(zed_muonsbis)")
            .Define("ssselected_muons",  "ReconstructedParticle::remove(sselected_muons, zed_muonsbis)")
            .Define("zed_muonsbisbis",   "ReconstructedParticle::findZleptons(ssselected_muons)")
            .Define("zed_muonic3",       "ReconstructedParticle::resonanceBuilder(91)(zed_muonsbisbis)")

#on fait une liste avec les muons qui n'ont pas pu former de paires, dans un cas de nombre impair 

            .Define("extramuons",        "ReconstructedParticle::remove(ssselected_muons, zed_muonsbisbis)")
            .Define("N_extramuons",      "ReconstructedParticle::get_n(extramuons)")

#même procédure que pour les muons mais avec les électrons

            .Define("zed_electrons",     "ReconstructedParticle::findZleptons(selected_electrons)")
            .Define("zed_electronic",    "ReconstructedParticle::resonanceBuilder(91)(zed_electrons)")
            .Define("sselected_electrons", "ReconstructedParticle::remove(selected_electrons, zed_electrons)")
            .Define("zed_electronsbis",      "ReconstructedParticle::findZleptons(sselected_electrons)")
            .Define("zed_electronic2",     "ReconstructedParticle::resonanceBuilder(91)(zed_electronsbis)")
            .Define("ssselected_electrons","ReconstructedParticle::remove(sselected_electrons, zed_electronsbis)")
            .Define("zed_electronsbisbis", "ReconstructedParticle::findZleptons(ssselected_electrons)")
            .Define("zed_electronic3",     "ReconstructedParticle::resonanceBuilder(91)(zed_electronsbisbis)")

#liste des électrons isolés/dans le range de p mais qui n'ont pas formé de paires

            .Define("extraelectrons",      "ReconstructedParticle::remove(ssselected_electrons, zed_electronsbisbis)")
            .Define("N_extraelectrons",    "ReconstructedParticle::get_n(extraelectrons)")

#on veut regrouper les 3 Z ensemble, donc les 3 muonic d'abord

            .Define("mergemuonic1",         "ReconstructedParticle::merge(zed_muonic, zed_muonic2)")
            .Define("mergemuonic2",         "ReconstructedParticle::merge(mergemuonic1, zed_muonic3)")

#et les 3 electronic (mais il n'y en aura qu'au maximum 3/6 non vides puisque trois Z max, ou alors bckg)
            
            .Define("mergeelectronic1",     "ReconstructedParticle::merge(zed_electronic, zed_electronic2)")
            .Define("mergeelectronic2",     "ReconstructedParticle::merge(mergeelectronic1, zed_electronic3)")

#et là on merge les 6 pour donner zed_leptonic 

            .Define("zed_leptonic",         "ReconstructedParticle::merge(mergemuonic2, mergeelectronic2)")

#on fait la liste de tous les leptons qui sont sélectionnés mais qui n'ont pas pu former de paires compatibles au Z

            .Define("extraleptons",         "ReconstructedParticle::merge(extramuons, extraelectrons)")

#puis on regarde le nombre de paires de Z et d'extraleptons

            .Define("N_zed_leptonic",       "ReconstructedParticle::get_n(zed_leptonic)")
            .Define("N_extraleptons",       "ReconstructedParticle::get_n(extraleptons)")
            

#propriétés de nos paires Z


            #get_e ne fonctionne pas sur zed_leptonic donc on peut le faire sur les leptons directement puis on sommera [0] et [1] pour avoir la composante d'énergie du tlv du Z (on n'a qu'un seul candidat au Z en stage 2 avec le fitlre donc on fait ça qu'avec les premières paires (et meilleures) d'électrons et muons de findZleptons 
            .Define("zed_electrons_e", "ReconstructedParticle::get_e(zed_electrons)")
            .Define("zed_muons_e", "ReconstructedParticle::get_e(zed_muons)")
            .Define("zed_leptons", "ReconstructedParticle::merge(zed_electrons, zed_muons)")
            .Define("zed_leptons_e", "ReconstructedParticle::get_e(zed_leptons)")
            
            .Define("zed_leptonic_m",       "ReconstructedParticle::get_mass(zed_leptonic)")
            .Define("zed_leptonic_pt",      "ReconstructedParticle::get_pt(zed_leptonic)")
            .Define("zed_leptonic_px",      "ReconstructedParticle::get_px(zed_leptonic)")
            .Define("zed_leptonic_py",      "ReconstructedParticle::get_py(zed_leptonic)") 
            .Define("zed_leptonic_pz",      "ReconstructedParticle::get_pz(zed_leptonic)")
            .Define("zed_leptonic_p",       "ReconstructedParticle::get_p(zed_leptonic)")
            .Define("zed_leptonic_recoil",  "ReconstructedParticle::recoilBuilder(240)(zed_leptonic)")
            .Define("zed_leptonic_recoil_m","ReconstructedParticle::get_mass(zed_leptonic_recoil)")            
            .Define("zed_leptonic_charge",  "ReconstructedParticle::get_charge(zed_leptonic)")
            .Define("zed_leptonic_theta",   "ReconstructedParticle::get_theta(zed_leptonic)")
            .Define("zed_leptonic_phi",     "ReconstructedParticle::get_phi(zed_leptonic)")
            .Define("zed_leptonic_y",       "ReconstructedParticle::get_y(zed_leptonic)")
            .Define("zed_leptonic_eta",     "ReconstructedParticle::get_eta(zed_leptonic)") 
            .Define("zed_leptonic_cos",     "cos(ReconstructedParticle::get_theta(zed_leptonic))")

#on fait la liste des leptons qui ont servi pour les paires : donc les selected - les extra
            
            .Define("taken_leptons",    "ReconstructedParticle::remove(selected_leptons, extraleptons)")
            .Define("N_taken_leptons",  "ReconstructedParticle::get_n(taken_leptons)")

#on sélectionne toutes les particules - celles utilisées pour les paires de Z pour reconstruire les jets
#puis on crée des pseudo-jets avec 

            .Define("my_recoparticles",  "ReconstructedParticle::remove(ReconstructedParticles, taken_leptons)")
            .Define("RP_px", "ReconstructedParticle::get_px(my_recoparticles)")
            .Define("RP_py", "ReconstructedParticle::get_py(my_recoparticles)")
            .Define("RP_pz", "ReconstructedParticle::get_pz(my_recoparticles)")
            .Define("RP_e", "ReconstructedParticle::get_e(my_recoparticles)")
            .Define("pseudo_jets",  "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")


            #Durham N=2

#construction des jets selon l'algo durham kt pour N=2 + (ligne2) mise en objet "pseudo_jets", dont on peut extraire des informations cinématiques
            .Define("FCCAnalysesJets_ee_genkt2",  "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
            .Define("jets_ee_genkt2",  "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt2)")
        



            .Define("jets_px2",  "JetClusteringUtils::get_px(jets_ee_genkt2)")
            .Define("jets_py2",  "JetClusteringUtils::get_py(jets_ee_genkt2)")
            .Define("jets_pz2",  "JetClusteringUtils::get_pz(jets_ee_genkt2)")
            .Define("jets_p2",   "JetClusteringUtils::get_p(jets_ee_genkt2)")
            .Define("jets_e2",   "JetClusteringUtils::get_e(jets_ee_genkt2)")
            .Define("jets_m2",   "JetClusteringUtils::get_m(jets_ee_genkt2)")
            .Define("jets_pt2",  "JetClusteringUtils::get_pt(jets_ee_genkt2)")
            .Define("jets_y2",   "JetClusteringUtils::get_y(jets_ee_genkt2)")
            .Define("jets_eta2", "JetClusteringUtils::get_eta(jets_ee_genkt2)")
            .Define("jets_theta2", "JetClusteringUtils::get_theta(jets_ee_genkt2)")
            .Define("jets_phi2", "JetClusteringUtils::get_phi(jets_ee_genkt2)")

            .Define("jetconstituents_ee_genkt2", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt2)")
            .Define("jetconstituents_ee_2", "JetConstituentsUtils::build_constituents_cluster(my_recoparticles, jetconstituents_ee_genkt2)")
            .Define("jetconstituents_2", "JetConstituentsUtils::count_consts(jetconstituents_ee_2)")
            
            

#j'ai l'impression que ça (la fonction get_constituents) donne le nombre de particules dans le jet 

            .Define("dmerge_2_45", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt2, 4)")
            .Define("dmerge_2_34", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt2, 3)")
            .Define("dmerge_2_23", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt2, 2)")
            .Define("dmerge_2_12", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt2, 1)")

            
            #.Define("2hadronic", "JetClusteringUtils::resonanceBuilder(91)(2jets_ee_genkt)")
            #.Define("2hadronic_mass", "JetClusteringUtils::get_m(2hadronic)")
            #.Define("2hadronic_e", "JetClusteringUtils::get_e(2hadronic)")
            #.Define("2hadronic_p", "JetClusteringUtils::get_p(2hadronic)")
            #.Define("2hadronic_pt", "JetClusteringUtils::get_pt(2hadronic)")


#"vraie" saveur reconstruite des jets (attention à la définition de get_flavour) (attention la VRAIE saveur c'est la saveur de la particule qui a donné naissance au jet qu'on connait via pdgid MONTE CARLO ; quand on reconstruit des jets on ne peut pas savoir laquelle exactement est à l'origine du jet donc on doit la déterminer au mieux via get_flavour qui n'est pas du tagging mais juste une manière d'attribuer une vraie valeur de saveur mais à un JET et non une particule, ce qui n'est pas évident)

            .Define("jets_ee_genkt_flavour2",   "JetTaggingUtils::get_flavour(jets_ee_genkt2, Particle)") #, 2, 0.8)")

            #nouveau get flavour : get_flavour_gm
            
#             .Define("jets_ee_flavour2", "JetTaggingUtils::get_flavour_gm(jets_ee_genkt2, Particle)")
            .Define("jets_ee_flavour2", "JetTaggingUtils::get_flavour(jets_ee_genkt2, Particle)")

#tagging (méthode pour les données obtenues par un détecteur) suivant la méthode paramétrée (voir ipad) : on définit des cuts nous-même pour l'algorithme, donnés avant le df dans ce fichier

            .Define("jets_ee_genkt_btag2",      "JetTaggingUtils::get_btag(jets_ee_genkt_flavour2, {:f}, {:f}, {:f}, {:f})".format(beff_b, beff_c, beff_l, beff_g))
            .Define("jets_ee_genkt_ctag2",      "JetTaggingUtils::get_ctag(jets_ee_genkt_flavour2, {:f}, {:f}, {:f}, {:f})".format(ceff_c, ceff_b, ceff_l, ceff_g))
            .Define("jets_ee_genkt_gtag2",      "JetTaggingUtils::get_gtag(jets_ee_genkt_flavour2, {:f}, {:f}, {:f}, {:f})".format(geff_g, geff_b, geff_c, geff_l))
            #.Define("jets_ee_genkt_ltag2",      "JetTaggingUtils::get_ltag(jets_ee_genkt_flavour2, {:f}, {:f}, {:f}, {:f})".format(leff_l, leff_b, leff_c, leff_g))

#nombre de jets taggés b,c,ou g suivant la méthode paramétrée ('fausse')

            .Define("n_jets_ee_genkt_btag2", "ReconstructedParticle::getJet_ntags(jets_ee_genkt_btag2)")
            .Define("n_jets_ee_genkt_ctag2", "ReconstructedParticle::getJet_ntags(jets_ee_genkt_ctag2)")
            .Define("n_jets_ee_genkt_gtag2", "ReconstructedParticle::getJet_ntags(jets_ee_genkt_gtag2)")
            #.Define("n_jets_ee_genkt_ltag2", "ReconstructedParticle::getJet_ntags(jets_ee_genkt_ltag2)")

#on renomme les pseudo jets parce que noms pénibles je crois

            .Alias("jets2", "jets_ee_genkt2")
            .Alias("jet_btag2", "jets_ee_genkt_btag2")
            .Alias("jet_ctag2", "jets_ee_genkt_ctag2")
            .Alias("jet_gtag2", "jets_ee_genkt_gtag2")
            #.Alias("jet_ltag2", "jets_ee_genkt_ltag2")

#ici on fait des listes de jets b,c,g et non b non c non g (avant on a taggé les jets mais ils sont pas groupés dans une liste)

            .Define("btagged_jets2", "JetTaggingUtils::sel_tag(true)(jet_btag2,jets2)")
            .Define("ctagged_jets2", "JetTaggingUtils::sel_tag(true)(jet_ctag2,jets2)")
            .Define("gtagged_jets2", "JetTaggingUtils::sel_tag(true)(jet_gtag2,jets2)")
            #.Define("ltagged_jets2", "JetTaggingUtils::sel_tag(true)(jet_ltag2,jets2)")

            .Define("nonbtagged_jets2", "JetTaggingUtils::sel_tag(false)(jet_btag2,jets2)")
            .Define("nonctagged_jets2", "JetTaggingUtils::sel_tag(false)(jet_ctag2,jets2)")
            .Define("nongtagged_jets2", "JetTaggingUtils::sel_tag(false)(jet_gtag2,jets2)")
            #.Define("nonltagged_jets2", "JetTaggingUtils::sel_tag(false)(jet_ltag2,jets2)")


#selection sur les jets
                
            .Define("selected_btagged_jets2",     "JetClusteringUtils::sel_p(0,150)(btagged_jets2)")
            .Define("selected_nonbtagged_jets2",  "JetClusteringUtils::sel_p(0,150)(nonbtagged_jets2)")
            .Define("selected_ctagged_jets2",     "JetClusteringUtils::sel_p(0,150)(ctagged_jets2)")
            .Define("selected_nonctagged_jets2",  "JetClusteringUtils::sel_p(0,150)(nonctagged_jets2)")
            .Define("selected_gtagged_jets2",     "JetClusteringUtils::sel_p(0,150)(gtagged_jets2)")
            .Define("selected_nongtagged_jets2",  "JetClusteringUtils::sel_p(0,150)(nongtagged_jets2)")
            #.Define("selected_ltagged_jets2",     "JetClusteringUtils::sel_p(0,100)(ltagged_jets2)")
            #.Define("selected_nonltagged_jets2", "JetClusteringUtils::sel_p(0,100)(nonltagged_jets2)")



#nombre de jets taggés b c ou g AVEC la sélection sur le momentum des jets

            .Define("n_selected_btagged_jets2",   "JetClusteringUtils::get_n(selected_btagged_jets2)")
            .Define("n_selected_ctagged_jets2",   "JetClusteringUtils::get_n(selected_ctagged_jets2)")
            .Define("n_selected_gtagged_jets2",   "JetClusteringUtils::get_n(selected_gtagged_jets2)")
            #.Define("n_selected_ltagged_jets2",  "JetClusteringUtils::get_n(selected_ltagged_jets2)")

#pt des jets sélectionnés en fonction de la saveur

            .Define("selected_btagged_jets_pt2",   "JetClusteringUtils::get_pt(selected_btagged_jets2)") 
            .Define("selected_nonbtagged_jets_pt2","JetClusteringUtils::get_pt(selected_nonbtagged_jets2)")
            .Define("selected_ctagged_jets_pt2",   "JetClusteringUtils::get_pt(selected_ctagged_jets2)") 
            .Define("selected_nonctagged_jets_pt2","JetClusteringUtils::get_pt(selected_nonctagged_jets2)")
            .Define("selected_gtagged_jets_pt2",   "JetClusteringUtils::get_pt(selected_gtagged_jets2)") 
            .Define("selected_nongtagged_jets_pt2","JetClusteringUtils::get_pt(selected_nongtagged_jets2)")
            #.Define("selected_ltagged_jets_pt2",   "JetClusteringUtils::get_pt(selected_ltagged_jets2)") 
            #.Define("selected_nonltagged_jets_pt2","JetClusteringUtils::get_pt(selected_nonltagged_jets2)")

#pareil mais pour p (normalement on devrait donc observer le cut)

            .Define("selected_btagged_jets_p2",    "JetClusteringUtils::get_p(selected_btagged_jets2)")
            .Define("selected_nonbtagged_jets_p2", "JetClusteringUtils::get_p(selected_nonbtagged_jets2)")
            .Define("selected_ctagged_jets_p2",    "JetClusteringUtils::get_p(selected_ctagged_jets2)")
            .Define("selected_nonctagged_jets_p2", "JetClusteringUtils::get_p(selected_nonctagged_jets2)")
            .Define("selected_gtagged_jets_p2",    "JetClusteringUtils::get_p(selected_gtagged_jets2)")
            .Define("selected_nongtagged_jets_p2", "JetClusteringUtils::get_p(selected_nongtagged_jets2)")
            #.Define("selected_ltagged_jets_p2",    "JetClusteringUtils::get_p(selected_ltagged_jets2)")
            #.Define("selected_nonltagged_jets_p2", "JetClusteringUtils::get_p(selected_nonltagged_jets2)")


#px

            .Define("selected_btagged_jets_px2",    "JetClusteringUtils::get_px(selected_btagged_jets2)")
            .Define("selected_nonbtagged_jets_px2", "JetClusteringUtils::get_px(selected_nonbtagged_jets2)")
            .Define("selected_ctagged_jets_px2",    "JetClusteringUtils::get_px(selected_ctagged_jets2)")
            .Define("selected_nonctagged_jets_px2", "JetClusteringUtils::get_px(selected_nonctagged_jets2)")
            .Define("selected_gtagged_jets_px2",    "JetClusteringUtils::get_px(selected_gtagged_jets2)")
            .Define("selected_nongtagged_jets_px2", "JetClusteringUtils::get_px(selected_nongtagged_jets2)")
            #.Define("selected_ltagged_jets_px2",    "JetClusteringUtils::get_px(selected_ltagged_jets2)")
            #.Define("selected_nonltagged_jets_px2", "JetClusteringUtils::get_px(selected_nonltagged_jets2)")


#py

            .Define("selected_btagged_jets_py2",    "JetClusteringUtils::get_py(selected_btagged_jets2)")
            .Define("selected_nonbtagged_jets_py2", "JetClusteringUtils::get_py(selected_nonbtagged_jets2)")
            .Define("selected_ctagged_jets_py2",    "JetClusteringUtils::get_py(selected_ctagged_jets2)")
            .Define("selected_nonctagged_jets_py2", "JetClusteringUtils::get_py(selected_nonctagged_jets2)")
            .Define("selected_gtagged_jets_py2",    "JetClusteringUtils::get_py(selected_gtagged_jets2)")
            .Define("selected_nongtagged_jets_py2", "JetClusteringUtils::get_py(selected_nongtagged_jets2)")

#pz

            .Define("selected_btagged_jets_pz2",    "JetClusteringUtils::get_pz(selected_btagged_jets2)")
            .Define("selected_nonbtagged_jets_pz2", "JetClusteringUtils::get_pz(selected_nonbtagged_jets2)")
            .Define("selected_ctagged_jets_pz2",    "JetClusteringUtils::get_pz(selected_ctagged_jets2)")
            .Define("selected_nonctagged_jets_pz2", "JetClusteringUtils::get_pz(selected_nonctagged_jets2)")
            .Define("selected_gtagged_jets_pz2",    "JetClusteringUtils::get_pz(selected_gtagged_jets2)")
            .Define("selected_nongtagged_jets_pz2", "JetClusteringUtils::get_pz(selected_nongtagged_jets2)")
            #.Define("selected_ltagged_jets_px2",    "JetClusteringUtils::get_px(selected_ltagged_jets2)")
            #.Define("selected_nonltagged_jets_px2", "JetClusteringUtils::get_px(selected_nonltagged_jets2)")


#pareil pour l'énergie

            .Define("selected_btagged_jets_e2",    "JetClusteringUtils::get_e(selected_btagged_jets2)")
            .Define("selected_nonbtagged_jets_e2", "JetClusteringUtils::get_e(selected_nonbtagged_jets2)")
            .Define("selected_ctagged_jets_e2",    "JetClusteringUtils::get_e(selected_ctagged_jets2)")
            .Define("selected_nonctagged_jets_e2", "JetClusteringUtils::get_e(selected_nonctagged_jets2)")
            .Define("selected_gtagged_jets_e2",    "JetClusteringUtils::get_e(selected_gtagged_jets2)")
            .Define("selected_nongtagged_jets_e2", "JetClusteringUtils::get_e(selected_nongtagged_jets2)")
            #.Define("selected_ltagged_jets_e2",    "JetClusteringUtils::get_e(selected_ltagged_jets2)")
            #.Define("selected_nonltagged_jets_e2", "JetClusteringUtils::get_e(selected_nonltagged_jets2)")



#pareil pour la masse

            .Define("selected_btagged_jets_m2",    "JetClusteringUtils::get_m(selected_btagged_jets2)")
            .Define("selected_nonbtagged_jets_m2", "JetClusteringUtils::get_m(selected_nonbtagged_jets2)")
            .Define("selected_ctagged_jets_m2",    "JetClusteringUtils::get_m(selected_ctagged_jets2)")
            .Define("selected_nonctagged_jets_m2", "JetClusteringUtils::get_m(selected_nonctagged_jets2)")
            .Define("selected_gtagged_jets_m2",    "JetClusteringUtils::get_m(selected_gtagged_jets2)")
            .Define("selected_nongtagged_jets_m2", "JetClusteringUtils::get_m(selected_nongtagged_jets2)")
            #.Define("selected_ltagged_jets_m2",    "JetClusteringUtils::get_m(selected_ltagged_jets2)")
            #.Define("selected_nonltagged_jets_m2", "JetClusteringUtils::get_m(selected_nonltagged_jets2)")
            




            
            
        


            
            
            
            

            #Durham N=4

            
            .Define("FCCAnalysesJets_ee_genkt4",  "JetClustering::clustering_ee_kt(2, 4, 1, 0)(pseudo_jets)")
            .Define("jets_ee_genkt4",  "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_genkt4)")
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

            .Define("jetconstituents_ee_genkt4", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_genkt4)")
            .Define("jetconstituents_ee_4", "JetConstituentsUtils::build_constituents_cluster(my_recoparticles, jetconstituents_ee_genkt4)"
)
            .Define("jetconstituents_4", "JetConstituentsUtils::count_consts(jetconstituents_ee_4)")
            .Define("jetconstituents_4_theta", "JetConstituentsUtils::get_theta(jetconstituents_ee_4)")
            .Define("jetconstituents_4_phi", "JetConstituentsUtils::get_phi(jetconstituents_ee_4)")
            .Define("jetconstituents_4_energy", "JetConstituentsUtils::get_e(jetconstituents_ee_4)")


            
            .Define("dmerge_4_45", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt4, 4)")
            .Define("dmerge_4_34", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt4, 3)")
            .Define("dmerge_4_23", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt4, 2)")
            .Define("dmerge_4_12", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_genkt4, 1)")




            #anti-kt R=0.4
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



            #fin anti-kt 

            #anti-kt R=0.6
            .Define("FCCAnalysesJets_ee_antikt4_bis", "JetClustering::clustering_antikt(0.6,0,0,0,0)(pseudo_jets)")
            .Define("jets_antikt4_bis", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_antikt4_bis)")
            .Define("N_jets_antikt4_bis", "JetClusteringUtils::get_n(jets_antikt4_bis)")
            .Define("jets_antikt_e4_bis", "JetClusteringUtils::get_e(jets_antikt4_bis)")
            .Define("jets_antikt_px4_bis", "JetClusteringUtils::get_px(jets_antikt4_bis)")
            .Define("jets_antikt_py4_bis", "JetClusteringUtils::get_py(jets_antikt4_bis)")
            .Define("jets_antikt_pz4_bis", "JetClusteringUtils::get_pz(jets_antikt4_bis)")
            .Define("jets_antikt_theta4_bis", "JetClusteringUtils::get_theta(jets_antikt4_bis)")
            .Define("jets_antikt_phi4_bis", "JetClusteringUtils::get_phi(jets_antikt4_bis)")

            .Define("jetconstituents_ee_antikt4_bis", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_antikt4_bis)")
            .Define("jetconstituents_ee_antikt4_utile_bis", "JetConstituentsUtils::build_constituents_cluster(my_recoparticles, jetconstituents_ee_antikt4_bis)")
            .Define("jetconstituents_antikt4_bis", "JetConstituentsUtils::count_consts(jetconstituents_ee_antikt4_utile_bis)")
            .Define("jetconstituents_antikt4_theta_bis", "JetConstituentsUtils::get_theta(jetconstituents_ee_antikt4_utile_bis)")
            .Define("jetconstituents_antikt4_phi_bis", "JetConstituentsUtils::get_phi(jetconstituents_ee_antikt4_utile_bis)")
            .Define("jetconstituents_antikt4_energy_bis", "JetConstituentsUtils::get_e(jetconstituents_ee_antikt4_utile_bis)")

            .Define("dmerge_antikt4_45_bis", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4_bis, 4)")
            .Define("dmerge_antikt4_34_bis", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4_bis, 3)")
            .Define("dmerge_antikt4_23_bis", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4_bis, 2)")
            .Define("dmerge_antikt4_12_bis", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4_bis, 1)")



            #fin anti-kt


             #anti-kt R=0.8
            .Define("FCCAnalysesJets_ee_antikt4_ter", "JetClustering::clustering_antikt(0.8,0,0,0,0)(pseudo_jets)")
            .Define("jets_antikt4_ter", "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_ee_antikt4_ter)")
            .Define("N_jets_antikt4_ter", "JetClusteringUtils::get_n(jets_antikt4_ter)")
            .Define("jets_antikt_e4_ter", "JetClusteringUtils::get_e(jets_antikt4_ter)")
            .Define("jets_antikt_px4_ter", "JetClusteringUtils::get_px(jets_antikt4_ter)")
            .Define("jets_antikt_py4_ter", "JetClusteringUtils::get_py(jets_antikt4_ter)")
            .Define("jets_antikt_pz4_ter", "JetClusteringUtils::get_pz(jets_antikt4_ter)")
            .Define("jets_antikt_theta4_ter", "JetClusteringUtils::get_theta(jets_antikt4_ter)")
            .Define("jets_antikt_phi4_ter", "JetClusteringUtils::get_phi(jets_antikt4_ter)")

            .Define("jetconstituents_ee_antikt4_ter", "JetClusteringUtils::get_constituents(FCCAnalysesJets_ee_antikt4_ter)")
            .Define("jetconstituents_ee_antikt4_utile_ter", "JetConstituentsUtils::build_constituents_cluster(my_recoparticles, jetconstituents_ee_antikt4_ter)")
            .Define("jetconstituents_antikt4_ter", "JetConstituentsUtils::count_consts(jetconstituents_ee_antikt4_utile_ter)")
            .Define("jetconstituents_antikt4_theta_ter", "JetConstituentsUtils::get_theta(jetconstituents_ee_antikt4_utile_ter)")
            .Define("jetconstituents_antikt4_phi_ter", "JetConstituentsUtils::get_phi(jetconstituents_ee_antikt4_utile_ter)")
            .Define("jetconstituents_antikt4_energy_ter", "JetConstituentsUtils::get_e(jetconstituents_ee_antikt4_utile_ter)")

            .Define("dmerge_antikt4_45_ter", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4_ter, 4)")
            .Define("dmerge_antikt4_34_ter", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4_ter, 3)")
            .Define("dmerge_antikt4_23_ter", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4_ter, 2)")
            .Define("dmerge_antikt4_12_ter", "JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_ee_antikt4_ter, 1)")



            #fin anti-kt

            
            #.Define("4hadronic", "ReconstructedParticle::multijetResonanceBuilder(0,2,3)(4jets_ee_genkt)")
            #.Define("4hadronic_mass", "JetClusteringUtils::get_m(4hadronic)")
            #.Define("4hadronic_e", "JetClusteringUtils::get_e(4hadronic)")
            #.Define("4hadronic_p", "JetClusteringUtils::get_p(4hadronic)")
            #.Define("4hadronic_pt", "JetClusteringUtils::get_pt(4hadronic)")

            # get truth jet flavour
        
#             .Define("jets_ee_genkt_flavour4",   "JetTaggingUtils::get_flavour(jets_ee_genkt4, Particle, 2, 0.8)")
            .Define("jets_ee_genkt_flavour4",   "JetTaggingUtils::get_flavour(jets_ee_genkt4, Particle)")   #, 2, 0.8)")


#             .Define("jets_ee_flavour4", "JetTaggingUtils::get_flavour_gm(jets_ee_genkt4, Particle)")
            .Define("jets_ee_flavour4", "JetTaggingUtils::get_flavour(jets_ee_genkt4, Particle)")
             
#check

             .Define("jets_ee_genkt_btag4",      "JetTaggingUtils::get_btag(jets_ee_genkt_flavour4, {:f}, {:f}, {:f}, {:f})".format(beff_b, beff_c, beff_l, beff_g))
            .Define("jets_ee_genkt_ctag4",      "JetTaggingUtils::get_ctag(jets_ee_genkt_flavour4, {:f}, {:f}, {:f}, {:f})".format(ceff_c, ceff_b, ceff_l, ceff_g))
            .Define("jets_ee_genkt_gtag4",      "JetTaggingUtils::get_gtag(jets_ee_genkt_flavour4, {:f}, {:f}, {:f}, {:f})".format(geff_g, geff_b, geff_c, geff_l))
            #.Define("jets_ee_genkt_ltag4",      "JetTaggingUtils::get_ltag(jets_ee_genkt_flavour4, {:f}, {:f}, {:f}, {:f})".format(leff_l, leff_b, leff_c, leff_g))



            .Define("n_jets_ee_genkt_btag4", "ReconstructedParticle::getJet_ntags(jets_ee_genkt_btag4)")
            .Define("n_jets_ee_genkt_ctag4", "ReconstructedParticle::getJet_ntags(jets_ee_genkt_ctag4)")
            .Define("n_jets_ee_genkt_gtag4", "ReconstructedParticle::getJet_ntags(jets_ee_genkt_gtag4)")
            #.Define("n_jets_ee_genkt_gtag4", "ReconstructedParticle::getJet_ntags(jets_ee_genkt_ltag4)")

            .Alias("jets4", "jets_ee_genkt4")
            .Alias("jet_btag4", "jets_ee_genkt_btag4")
            .Alias("jet_ctag4", "jets_ee_genkt_ctag4")
            .Alias("jet_gtag4", "jets_ee_genkt_gtag4")
            #.Alias("jet_ltag4", "jets_ee_genkt_ltag4")

            .Define("btagged_jets4", "JetTaggingUtils::sel_tag(true)(jet_btag4,jets4)")
            .Define("ctagged_jets4", "JetTaggingUtils::sel_tag(true)(jet_ctag4,jets4)")
            .Define("gtagged_jets4", "JetTaggingUtils::sel_tag(true)(jet_gtag4,jets4)")
            #.Define("ltagged_jets4", "JetTaggingUtils::sel_tag(true)(jet_ltag4,jets4)")
            .Define("nonbtagged_jets4", "JetTaggingUtils::sel_tag(false)(jet_btag4,jets4)")
            .Define("nonctagged_jets4", "JetTaggingUtils::sel_tag(false)(jet_ctag4,jets4)")
            .Define("nongtagged_jets4", "JetTaggingUtils::sel_tag(false)(jet_gtag4,jets4)")
            #.Define("nonltagged_jets4", "JetTaggingUtils::sel_tag(false)(jets_ltag4,jets4)")
            
                
            .Define("selected_btagged_jets4",     "JetClusteringUtils::sel_p(0,100)(btagged_jets4)")
            .Define("selected_nonbtagged_jets4",  "JetClusteringUtils::sel_p(0,100)(nonbtagged_jets4)")
            .Define("selected_ctagged_jets4",     "JetClusteringUtils::sel_p(0,100)(ctagged_jets4)")
            .Define("selected_nonctagged_jets4",  "JetClusteringUtils::sel_p(0,100)(nonctagged_jets4)")
            .Define("selected_gtagged_jets4",     "JetClusteringUtils::sel_p(0,100)(gtagged_jets4)")
            .Define("selected_nongtagged_jets4",  "JetClusteringUtils::sel_p(0,100)(nongtagged_jets4)")
            #.Define("selected_ltagged_jets4",     "JetClusteringUtils::sel_p(0,100)(ltagged_jets4)")
            #.Define("selected_nonltagged_jets4",  "JetClusteringUtils::sel_p(0,100)(nonltagged_jets4)")


            .Define("n_selected_btagged_jets4",   "JetClusteringUtils::get_n(selected_btagged_jets4)")
            .Define("n_selected_ctagged_jets4",   "JetClusteringUtils::get_n(selected_ctagged_jets4)")
            .Define("n_selected_gtagged_jets4",   "JetClusteringUtils::get_n(selected_gtagged_jets4)")
            #.Define("n_selected_ltagged_jets4",   "JetClusteringUtils::get_n(selected_ltagged_jets4)")

            .Define("selected_btagged_jets_pt4",   "JetClusteringUtils::get_pt(selected_btagged_jets4)") 
            .Define("selected_nonbtagged_jets_pt4","JetClusteringUtils::get_pt(selected_nonbtagged_jets4)")
            .Define("selected_ctagged_jets_pt4",   "JetClusteringUtils::get_pt(selected_ctagged_jets4)") 
            .Define("selected_nonctagged_jets_pt4","JetClusteringUtils::get_pt(selected_nonctagged_jets4)")
            .Define("selected_gtagged_jets_pt4",   "JetClusteringUtils::get_pt(selected_gtagged_jets4)") 
            .Define("selected_nongtagged_jets_pt4","JetClusteringUtils::get_pt(selected_nongtagged_jets4)")
            #.Define("selected_ltagged_jets_pt4",   "JetClusteringUtils::get_pt(selected_ltagged_jets4)") 
            #.Define("selected_nonltagged_jets_pt4","JetClusteringUtils::get_pt(selected_nonltagged_jets4)")

            .Define("selected_btagged_jets_p4",    "JetClusteringUtils::get_p(selected_btagged_jets4)")
            .Define("selected_nonbtagged_jets_p4", "JetClusteringUtils::get_p(selected_nonbtagged_jets4)")
            .Define("selected_ctagged_jets_p4",    "JetClusteringUtils::get_p(selected_ctagged_jets4)")
            .Define("selected_nonctagged_jets_p4", "JetClusteringUtils::get_p(selected_nonctagged_jets4)")
            .Define("selected_gtagged_jets_p4",    "JetClusteringUtils::get_p(selected_gtagged_jets4)")
            .Define("selected_nongtagged_jets_p4", "JetClusteringUtils::get_p(selected_nongtagged_jets4)")
            #.Define("selected_ltagged_jets_p4",    "JetClusteringUtils::get_p(selected_ltagged_jets4)")
            #.Define("selected_nonltagged_jets_p4", "JetClusteringUtils::get_p(selected_nonltagged_jets4)")

            .Define("selected_btagged_jets_px4",    "JetClusteringUtils::get_px(selected_btagged_jets4)")
            .Define("selected_nonbtagged_jets_px4", "JetClusteringUtils::get_px(selected_nonbtagged_jets4)")
            .Define("selected_ctagged_jets_px4",    "JetClusteringUtils::get_px(selected_ctagged_jets4)")
            .Define("selected_nonctagged_jets_px4", "JetClusteringUtils::get_px(selected_nonctagged_jets4)")
            .Define("selected_gtagged_jets_px4",    "JetClusteringUtils::get_px(selected_gtagged_jets4)")
            .Define("selected_nongtagged_jets_px4", "JetClusteringUtils::get_px(selected_nongtagged_jets4)")
            #.Define("selected_ltagged_jets_px4",    "JetClusteringUtils::get_px(selected_ltagged_jets4)")
            #.Define("selected_nonltagged_jets_px4", "JetClusteringUtils::get_px(selected_nonltagged_jets4)")

            .Define("selected_btagged_jets_py4",    "JetClusteringUtils::get_py(selected_btagged_jets4)")
            .Define("selected_nonbtagged_jets_py4", "JetClusteringUtils::get_py(selected_nonbtagged_jets4)")
            .Define("selected_ctagged_jets_py4",    "JetClusteringUtils::get_py(selected_ctagged_jets4)")
            .Define("selected_nonctagged_jets_py4", "JetClusteringUtils::get_py(selected_nonctagged_jets4)")
            .Define("selected_gtagged_jets_py4",    "JetClusteringUtils::get_py(selected_gtagged_jets4)")
            .Define("selected_nongtagged_jets_py4", "JetClusteringUtils::get_py(selected_nongtagged_jets4)")
            #.Define("selected_ltagged_jets_py4",    "JetClusteringUtils::get_py(selected_ltagged_jets4)")
            #.Define("selected_nonltagged_jets_py4", "JetClusteringUtils::get_py(selected_nonltagged_jets4)")

            .Define("selected_btagged_jets_pz4",    "JetClusteringUtils::get_pz(selected_btagged_jets4)")
            .Define("selected_nonbtagged_jets_pz4", "JetClusteringUtils::get_pz(selected_nonbtagged_jets4)")
            .Define("selected_ctagged_jets_pz4",    "JetClusteringUtils::get_pz(selected_ctagged_jets4)")
            .Define("selected_nonctagged_jets_pz4", "JetClusteringUtils::get_pz(selected_nonctagged_jets4)")
            .Define("selected_gtagged_jets_pz4",    "JetClusteringUtils::get_pz(selected_gtagged_jets4)")
            .Define("selected_nongtagged_jets_pz4", "JetClusteringUtils::get_pz(selected_nongtagged_jets4)")
            #.Define("selected_ltagged_jets_pz4",    "JetClusteringUtils::get_pz(selected_ltagged_jets4)")
            #.Define("selected_nonltagged_jets_pz4", "JetClusteringUtils::get_pz(selected_nonltagged_jets4)")
            

            .Define("selected_btagged_jets_e4",    "JetClusteringUtils::get_e(selected_btagged_jets4)")
            .Define("selected_nonbtagged_jets_e4", "JetClusteringUtils::get_e(selected_nonbtagged_jets4)")
            .Define("selected_ctagged_jets_e4",    "JetClusteringUtils::get_e(selected_ctagged_jets4)")
            .Define("selected_nonctagged_jets_e4", "JetClusteringUtils::get_e(selected_nonctagged_jets4)")
            .Define("selected_gtagged_jets_e4",    "JetClusteringUtils::get_e(selected_gtagged_jets4)")
            .Define("selected_nongtagged_jets_e4", "JetClusteringUtils::get_e(selected_nongtagged_jets4)")
            #.Define("selected_ltagged_jets_e4",    "JetClusteringUtils::get_e(selected_ltagged_jets4)")
            #.Define("selected_nonltagged_jets_e4", "JetClusteringUtils::get_e(selected_nonltagged_jets4)")

            .Define("selected_btagged_jets_m4",    "JetClusteringUtils::get_m(selected_btagged_jets4)")
            .Define("selected_nonbtagged_jets_m4", "JetClusteringUtils::get_m(selected_nonbtagged_jets4)")
            .Define("selected_ctagged_jets_m4",    "JetClusteringUtils::get_m(selected_ctagged_jets4)")
            .Define("selected_nonctagged_jets_m4", "JetClusteringUtils::get_m(selected_nonctagged_jets4)")
            .Define("selected_gtagged_jets_m4",    "JetClusteringUtils::get_m(selected_gtagged_jets4)")
            .Define("selected_nongtagged_jets_m4", "JetClusteringUtils::get_m(selected_nongtagged_jets4)")
            #.Define("selected_ltagged_jets_m4",    "JetClusteringUtils::get_m(selected_ltagged_jets4)")
            #.Define("selected_nonltagged_jets_m4", "JetClusteringUtils::get_m(selected_nonltagged_jets4)")

            

            
            
            

        

#hzz monte carlo

            .Alias("Particle1", "Particle#1.index")
#             .Define("hzz_decay", "MCParticle::fill_ZHZZ_decay(Particle, Particle1)")
            .Define("ZH_decay", "MCParticle::fill_ZH_decay(Particle, Particle1)")
#             .Define("inv_mass_H", "MCParticle::invariant_mass(Hbb_decay.Z_decay)")

#             .Define("inv_mass_Z", "MCParticle::invariant_mass(hzz_decay.Z_decay)")
#             .Define("pdg_Z", "MCParticle::get_pdg(hzz_decay.Z_decay)")

#             .Define("inv_mass_Z1", "MCParticle::invariant_mass(hzz_decay.Z1_decay)")
#             .Define("pdg_Z1", "MCParticle::get_pdg(hzz_decay.Z1_decay)")

#             .Define("inv_mass_Z2", "MCParticle::invariant_mass(hzz_decay.Z2_decay)")
#             .Define("pdg_Z2", "MCParticle::get_pdg(hzz_decay.Z2_decay)")

#testé ici mais va au stage 2 Zleptonique
            #.Define("Zelectronique", "MCParticle::sel_pdgID(11, true)(hzz_decay.Z_decay)")
            #.Define("Zmuonique","MCParticle::sel_pdgID(13, true)(hzz_decay.Z_decay)")
            #.Define("Zleptonique", "MCParticle::mergeParticles(Zelectronique, Zmuonique)")

            
#prendre les MCParticles qui ont passé le filtre de Z leptonique : liste Zleptonique 
#il faut faire pareil que pour Z mais pour Z1 et Z2 (l'un après l'autre c'est ok je pense à chaque fois faut juste prendre la liste qui a subi les sélections) , mais là attention il faut faire un sel_pdgID avec des sortes de où pour avoir tous les quarks/gluon : écrire une fonction ? en se basant sur sel_pdgID mais en rajoutant + d'arguments en entrée et un ou dans les condit

            .Define("emiss", "MissingET.energy[0]")
            .Define("pxmiss", "MissingET.momentum.x[0]")
            .Define("pymiss", "MissingET.momentum.y[0]")
            .Define("pzmiss", "MissingET.momentum.z[0]")
            .Define("missing_tlv", "ReconstructedParticle::get_tlv(MissingET)")
            .Define("etmiss", "sqrt((MissingET.momentum.x[0])*(MissingET.momentum.x[0]) + (MissingET.momentum.y[0])*(MissingET.momentum.y[0]))")

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

            "jets_px2",
            "jets_py2",
            "jets_pz2",
            "jets_p2",
            "jets_e2",
            "jets_m2",
            "jets_pt2",
            "jets_y2",
            "jets_eta2",
            "jets_theta2",
            "jets_phi2",
            "jetconstituents_ee_genkt2",
            "jetconstituents_ee_2",
            "jetconstituents_2",
            "dmerge_2_45",
            "dmerge_2_34",
            "dmerge_2_23",
            "dmerge_2_12",

            #"2hadronic_mass",
            #"2hadronic_e",
            #"2hadronic_p",
            #"2hadronic_pt",

            "jets_ee_genkt_flavour2",
            "jets_ee_flavour2",

            "n_jets_ee_genkt_btag2",
            "n_jets_ee_genkt_ctag2",
            "n_jets_ee_genkt_gtag2",
            "n_selected_btagged_jets2",
            "n_selected_ctagged_jets2",
            "n_selected_gtagged_jets2",

            "selected_btagged_jets_pt2",
            "selected_nonbtagged_jets_pt2",
            "selected_ctagged_jets_pt2",
            "selected_nonctagged_jets_pt2",
            "selected_gtagged_jets_pt2",
            "selected_nongtagged_jets_pt2",

            

            "selected_btagged_jets_pz2",
            "selected_nonbtagged_jets_pz2",
            "selected_ctagged_jets_pz2",
            "selected_nonctagged_jets_pz2",
            "selected_gtagged_jets_pz2",
            "selected_nongtagged_jets_pz2",

            "selected_btagged_jets_p2",
            "selected_nonbtagged_jets_p2",
            "selected_ctagged_jets_p2",
            "selected_nonctagged_jets_p2",
            "selected_gtagged_jets_p2",
            "selected_nongtagged_jets_p2",

            "selected_btagged_jets_e2",
            "selected_nonbtagged_jets_e2",
            "selected_ctagged_jets_e2",
            "selected_nonctagged_jets_e2",
            "selected_gtagged_jets_e2",
            "selected_nongtagged_jets_e2",

            "selected_btagged_jets_m2",
            "selected_nonbtagged_jets_m2",
            "selected_ctagged_jets_m2",
            "selected_nonctagged_jets_m2",
            "selected_gtagged_jets_m2",
            "selected_nongtagged_jets_m2",

            
    
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

            #"4hadronic_mass",
            #"4hadronic_e",
            #"4hadronic_p",
            #"4hadronic_pt",

            "jets_ee_genkt_flavour4",
            "jets_ee_flavour4",
            "jetconstituents_4_theta",
            "jetconstituents_4_phi",
            "jetconstituents_4_energy",

            
            "n_jets_ee_genkt_btag4",
            "n_jets_ee_genkt_ctag4",
            "n_jets_ee_genkt_gtag4",
            "n_selected_btagged_jets4",
            "n_selected_ctagged_jets4",
            "n_selected_gtagged_jets4",

            "selected_btagged_jets_pt4",
            "selected_nonbtagged_jets_pt4",
            "selected_ctagged_jets_pt4",
            "selected_nonctagged_jets_pt4",
            "selected_gtagged_jets_pt4",
            "selected_nongtagged_jets_pt4",

            "selected_btagged_jets_px4",
            "selected_nonbtagged_jets_px4",
            "selected_ctagged_jets_px4",
            "selected_nonctagged_jets_px4",
            "selected_gtagged_jets_px4",
            "selected_nongtagged_jets_px4",

            "selected_btagged_jets_py4",
            "selected_nonbtagged_jets_py4",
            "selected_ctagged_jets_py4",
            "selected_nonctagged_jets_py4",
            "selected_gtagged_jets_py4",
            "selected_nongtagged_jets_py4",

            "selected_btagged_jets_pz4",
            "selected_nonbtagged_jets_pz4",
            "selected_ctagged_jets_pz4",
            "selected_nonctagged_jets_pz4",
            "selected_gtagged_jets_pz4",
            "selected_nongtagged_jets_pz4",

            "selected_btagged_jets_p4",
            "selected_nonbtagged_jets_p4",
            "selected_ctagged_jets_p4",
            "selected_nonctagged_jets_p4",
            "selected_gtagged_jets_p4",
            "selected_nongtagged_jets_p4",

            "selected_btagged_jets_e4",
            "selected_nonbtagged_jets_e4",
            "selected_ctagged_jets_e4",
            "selected_nonctagged_jets_e4",
            "selected_gtagged_jets_e4",
            "selected_nongtagged_jets_e4",

            "selected_btagged_jets_m4",
            "selected_nonbtagged_jets_m4",
            "selected_ctagged_jets_m4",
            "selected_nonctagged_jets_m4",
            "selected_gtagged_jets_m4",
            "selected_nongtagged_jets_m4",

    


            
            
#             "hzz_decay",
            "ZH_decay",
#             "inv_mass_H",
#             "inv_mass_Z",
#             "inv_mass_H",
#             "pdg_Z",

#             "inv_mass_Z1",
#             "pdg_Z1",

#             "inv_mass_Z2",
#             "pdg_Z2",

            "emiss",
            "pxmiss",
            "pymiss",
            "pzmiss",
            "etmiss",
	    "missing_tlv",


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

            "N_jets_antikt4_bis",
            "jets_antikt_e4_bis",
            "jets_antikt_px4_bis",
            "jets_antikt_py4_bis",
            "jets_antikt_pz4_bis",
            "jets_antikt_theta4_bis",
            "jets_antikt_phi4_bis",
            "jetconstituents_antikt4_bis",
            "jetconstituents_antikt4_theta_bis",
            "jetconstituents_antikt4_phi_bis",
            "jetconstituents_antikt4_energy_bis",
            "dmerge_antikt4_45_bis",
            "dmerge_antikt4_34_bis",
            "dmerge_antikt4_23_bis",
            "dmerge_antikt4_12_bis",

            "N_jets_antikt4_ter",
            "jets_antikt_e4_ter",
            "jets_antikt_px4_ter",
            "jets_antikt_py4_ter",
            "jets_antikt_pz4_ter",
            "jets_antikt_theta4_ter",
            "jets_antikt_phi4_ter",
            "jetconstituents_antikt4_ter",
            "jetconstituents_antikt4_theta_ter",
            "jetconstituents_antikt4_phi_ter",
            "jetconstituents_antikt4_energy_ter",
            "dmerge_antikt4_45_ter",
            "dmerge_antikt4_34_ter",
            "dmerge_antikt4_23_ter",
            "dmerge_antikt4_12_ter",

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


