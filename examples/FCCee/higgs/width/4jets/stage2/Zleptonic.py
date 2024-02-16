inputDir    = "/usatlas/u/atishelma/FCC/FCCAnalyses/2DPlotting/stage1"

#Optional: output directory, default is local dir
outputDir    = "/usatlas/u/atishelma/FCC/FCCAnalyses/2DPlotting/stage2"
# outputDir   = "outputs/fccee/higgs/mH-recoil/hzz/stage2/4jets/Zleptonic/"

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
              # 'wzp6_ee_eeH_HZZ_ecm240':{}
               #'wzp6_ee_eeH_HWW_ecm240':{},
               #'wzp6_ee_eeH_HZa_ecm240':{},
              # 'wzp6_ee_eeH_Haa_ecm240':{},
              # 'wzp6_ee_eeH_Hbb_ecm240':{},
              # 'wzp6_ee_eeH_Hcc_ecm240':{},
              # 'wzp6_ee_eeH_Hgg_ecm240':{},
              # 'wzp6_ee_eeH_Hmumu_ecm240':{},
              # 'wzp6_ee_eeH_Hss_ecm240':{},
              # 'wzp6_ee_eeH_Htautau_ecm240':{},
              # 'wzp6_ee_ee_Mee_30_150_ecm240':{},
              # 'p8_ee_ZZ_ecm240':{},#Run the full statistics in one output file named <outputDir>/p8_ee_ZZ_ecm240.root
#                'p8_ee_WW_ecm240':{} #Run 50% of the statistics in two files named <outputDir>/p8_ee_WW_ecm240/chunk<N>.root
    #'p8_ee_ZH_ecm240':{'fraction':0.2, 'output':'p8_ee_ZH_ecm240_out'} #Run 20% of the statistics in one file named <outputDir>/p8_ee_ZH_ecm240_out.root (example on how to change the output name)
               }


#Optional: ncpus, default is 4
nCPUS       = 64 

#Optional running on HTCondor, default is False
#runBatch    = False

#USER DEFINED CODE
import ROOT
ROOT.gInterpreter.Declare("""
bool myFilter(ROOT::VecOps::RVec<float> mass) {
    for (size_t i = 0; i < mass.size(); ++i) {
        if (mass.at(i)>80. && mass.at(i)<100.)
            return true;
    }
    return false;
}
""")
#END USER DEFINED CODE

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (df
               #pour tourner que sur du signal donc sur du hzz
               #.Filter("hzz_decay.Z_decay.size()>0")

               #on stocke la valeur ABSOLUE du pdg de Z ; attention, pdg_Z stocke 2 variables donc il faut en choisir une avant d'appliquer le filtre, en l'occurence pour hzz truth les deux valeurs sont juste de signe opposé donc on prend la valeur absolue et voilà
               #.Define("abs_pdg_Z", "abs(MCParticle::get_pdg(hzz_decay.Z_decay)[0])")
               #on applique le filtre pour avoir un decay du Z seul en leptons (électrons ou muons)
              # .Filter("abs_pdg_Z == 11 || abs_pdg_Z == 13")
               
               #même chose mais maintenant pour les deux Z du Higgs qui eux sont voulus en jets (quarks/gluons en truth)

               #.Define("abs_pdg_Z1", "abs(MCParticle::get_pdg(hzz_decay.Z1_decay)[0])")
               #.Filter("abs_pdg_Z1 == 1 || abs_pdg_Z1 == 2 || abs_pdg_Z1 == 3 || abs_pdg_Z1 == 4 || abs_pdg_Z1 == 5")
               #.Define("abs_pdg_Z2", "abs(MCParticle::get_pdg(hzz_decay.Z2_decay)[0])")
               #.Filter("abs_pdg_Z2 == 1 || abs_pdg_Z2 == 2 || abs_pdg_Z2 == 3 || abs_pdg_Z2 == 4 || abs_pdg_Z2 == 5")

        
               
#maintenant, on est plus que sur des évènements où le Z seul est en leptons et le higgs en deux Z eux meme en 2 jets (4jets)

               .Filter("size(jets_e4)>0") 
               #Filter to have exactly one Z candidate for the reconstructed particles
#                .Filter("zed_leptonic_m.size() == 1")
               #.Define Z candidate mass
               .Define("Zcand_m","zed_leptonic_m[0]")
               #Define Z candidate recoil mass
               .Define("Zcand_recoil_m","zed_leptonic_recoil_m[0]")
               #Define Z candidate pt
               .Define("Zcand_pt","zed_leptonic_pt[0]")
               .Define("Zcand_px","zed_leptonic_px[0]")
               .Define("Zcand_py","zed_leptonic_py[0]")
               .Define("Zcand_pz","zed_leptonic_pz[0]")
               .Define("Zcand_p","zed_leptonic_p[0]")
               .Define("Zcand_q","zed_leptonic_charge[0]")
               .Define("Zcand_phi", "zed_leptonic_phi[0]")
               .Define("Zcand_theta", "zed_leptonic_theta[0]")
               .Define("Zcand_cos","zed_leptonic_cos[0]")
               .Define("Zcand_y", "zed_leptonic_y[0]")
               .Define("Zcand_eta", "zed_leptonic_eta[0]")

               #on considère les 4 jets de Duhram

               #resonance Builder: on regroupe les jets par 2 ! 

               #.Define("the2Z", "ReconstructedParticle::myresoBuilder(jets_e4, jets_px4, jets_py4, jets_pz4, jets_ee_genkt_flavour4)")
               .Define("the2Z", "ReconstructedParticle::myresoBuilder(jets_e4, jets_px4, jets_py4, jets_pz4, jets_ee_genkt_flavour4, jets_ee_flavour4, jets_eta4, jets_theta4, jets_phi4)")

               .Define("jetconstituents_4_1", "jetconstituents_4[0]")
               .Define("jetconstituents_4_2", "jetconstituents_4[1]")
               .Define("jetconstituents_4_3", "jetconstituents_4[2]")
               .Define("jetconstituents_4_4", "jetconstituents_4[3]")

               #duhram antikt R=0.4

               .Define("N_jets_antikt4_Ab5", "ReconstructedParticle::countNjets(jets_antikt_e4, 5)")
               .Define("N_jets_antikt4_Ab10", "ReconstructedParticle::countNjets(jets_antikt_e4, 10)")
               

               .Define("jetconstituents_antikt4_1", "jetconstituents_antikt4[0]")
               .Define("jetconstituents_antikt4_2", "jetconstituents_antikt4[1]")
               .Define("jetconstituents_antikt4_3", "jetconstituents_antikt4[2]")
               .Define("jetconstituents_antikt4_4", "jetconstituents_antikt4[3]")

               .Define("jet1_antikt4_energy", "jets_antikt_e4[0]")
               .Define("jet2_antikt4_energy", "jets_antikt_e4[1]")
               .Define("jet3_antikt4_energy", "jets_antikt_e4[2]")
               .Define("jet4_antikt4_energy", "jets_antikt_e4[3]")
               .Define("jet5_antikt4_energy", "jets_antikt_e4[4]")
               .Define("jet6_antikt4_energy", "jets_antikt_e4[5]")
               .Define("jet7_antikt4_energy", "jets_antikt_e4[6]")
               .Define("jet8_antikt4_energy", "jets_antikt_e4[7]")
               .Define("jet9_antikt4_energy", "jets_antikt_e4[8]")
               .Define("jet10_antikt4_energy", "jets_antikt_e4[9]")
               .Define("jet11_antikt4_energy", "jets_antikt_e4[10]")
               .Define("jet12_antikt4_energy", "jets_antikt_e4[11]")

               .Define("jetconstituents_antikt4_theta_1", "jetconstituents_antikt4_theta[0]")
               .Define("jetconstituents_antikt4_theta_2", "jetconstituents_antikt4_theta[1]")
               .Define("jetconstituents_antikt4_theta_3", "jetconstituents_antikt4_theta[2]")
               .Define("jetconstituents_antikt4_theta_4", "jetconstituents_antikt4_theta[3]")

               .Define("jetconstituents_antikt4_phi_1", "jetconstituents_antikt4_phi[0]")
               .Define("jetconstituents_antikt4_phi_2", "jetconstituents_antikt4_phi[1]")
               .Define("jetconstituents_antikt4_phi_3", "jetconstituents_antikt4_phi[2]")
               .Define("jetconstituents_antikt4_phi_4", "jetconstituents_antikt4_phi[3]")

               .Define("jetconstituents_antikt4_energy_1", "jetconstituents_antikt4_energy[0]")
               .Define("jetconstituents_antikt4_energy_2", "jetconstituents_antikt4_energy[1]")
               .Define("jetconstituents_antikt4_energy_3", "jetconstituents_antikt4_energy[2]")
               .Define("jetconstituents_antikt4_energy_4", "jetconstituents_antikt4_energy[3]")

               #Reconstruction des Z à partir des jets

               .Define("the2Z_antikt", "ReconstructedParticle::resoantikt(jets_antikt_e4, jets_antikt_px4, jets_antikt_py4 ,jets_antikt_pz4, N_jets_antikt4_Ab5, jets_antikt_theta4)")

               #methode : Z1 : 1+2+5+6+... ; Z2 = 3 + 4

               .Define("firstZ_antikt", "the2Z_antikt.Z1")
               .Define("firstZ_m_antikt" ,"firstZ_antikt.M()")

               .Define("secondZ_antikt", "the2Z_antikt.Z2")
               .Define("secondZ_m_antikt","secondZ_antikt.M()")

               #methode : Z1 = 1+2

               .Define("firstZ_antikt_1_2", "the2Z_antikt.Z1_1_2")
               .Define("firstZ_m_antikt_1_2", "firstZ_antikt_1_2.M()")

               #methode : Z1 = 1+2 + tout ce qui est + proche de 1 ou 2 que de 3 ou 4 ; Z2 = 3+4+ ce qui n'est pas avec Z1

               .Define("firstZ_antikt_reco", "the2Z_antikt.Z1_reco")
               .Define("firstZ_m_antikt_reco", "firstZ_antikt_reco.M()")

               #fin

               #duhram antikt R=0.6

               .Define("N_jets_antikt4_Ab5_bis", "ReconstructedParticle::countNjets(jets_antikt_e4_bis, 5)")
               .Define("N_jets_antikt4_Ab10_bis", "ReconstructedParticle::countNjets(jets_antikt_e4_bis, 10)")

               .Define("jetconstituents_antikt4_1_bis", "jetconstituents_antikt4_bis[0]")
               .Define("jetconstituents_antikt4_2_bis", "jetconstituents_antikt4_bis[1]")
               .Define("jetconstituents_antikt4_3_bis", "jetconstituents_antikt4_bis[2]")
               .Define("jetconstituents_antikt4_4_bis", "jetconstituents_antikt4_bis[3]")
               .Define("jetconstituents_antikt4_5_bis", "jetconstituents_antikt4_bis[4]")
               .Define("jetconstituents_antikt4_6_bis", "jetconstituents_antikt4_bis[5]")



               .Define("jet1_antikt4_energy_bis", "jets_antikt_e4_bis[0]")
               .Define("jet2_antikt4_energy_bis", "jets_antikt_e4_bis[1]")
               .Define("jet3_antikt4_energy_bis", "jets_antikt_e4_bis[2]")
               .Define("jet4_antikt4_energy_bis", "jets_antikt_e4_bis[3]")
               .Define("jet5_antikt4_energy_bis", "jets_antikt_e4_bis[4]")
               .Define("jet6_antikt4_energy_bis", "jets_antikt_e4_bis[5]")
               .Define("jet7_antikt4_energy_bis", "jets_antikt_e4_bis[6]")
               .Define("jet8_antikt4_energy_bis", "jets_antikt_e4_bis[7]")
               .Define("jet9_antikt4_energy_bis", "jets_antikt_e4_bis[8]")
               .Define("jet10_antikt4_energy_bis", "jets_antikt_e4_bis[9]")
               .Define("jet11_antikt4_energy_bis", "jets_antikt_e4_bis[10]")
               .Define("jet12_antikt4_energy_bis", "jets_antikt_e4_bis[11]")

               .Define("jetconstituents_antikt4_theta_1_bis", "jetconstituents_antikt4_theta_bis[0]")
               .Define("jetconstituents_antikt4_theta_2_bis", "jetconstituents_antikt4_theta_bis[1]")
               .Define("jetconstituents_antikt4_theta_3_bis", "jetconstituents_antikt4_theta_bis[2]")
               .Define("jetconstituents_antikt4_theta_4_bis", "jetconstituents_antikt4_theta_bis[3]")               
               .Define("jetconstituents_antikt4_theta_5_bis", "jetconstituents_antikt4_theta_bis[4]")
               .Define("jetconstituents_antikt4_theta_6_bis", "jetconstituents_antikt4_theta_bis[5]")


               .Define("jetconstituents_antikt4_phi_1_bis", "jetconstituents_antikt4_phi_bis[0]")
               .Define("jetconstituents_antikt4_phi_2_bis", "jetconstituents_antikt4_phi_bis[1]")
               .Define("jetconstituents_antikt4_phi_3_bis", "jetconstituents_antikt4_phi_bis[2]")
               .Define("jetconstituents_antikt4_phi_4_bis", "jetconstituents_antikt4_phi_bis[3]")
               .Define("jetconstituents_antikt4_phi_5_bis", "jetconstituents_antikt4_phi_bis[4]")
               .Define("jetconstituents_antikt4_phi_6_bis", "jetconstituents_antikt4_phi_bis[5]")

               .Define("jetconstituents_antikt4_energy_1_bis", "jetconstituents_antikt4_energy_bis[0]")
               .Define("jetconstituents_antikt4_energy_2_bis", "jetconstituents_antikt4_energy_bis[1]")
               .Define("jetconstituents_antikt4_energy_3_bis", "jetconstituents_antikt4_energy_bis[2]")
               .Define("jetconstituents_antikt4_energy_4_bis", "jetconstituents_antikt4_energy_bis[3]")
               .Define("jetconstituents_antikt4_energy_5_bis", "jetconstituents_antikt4_energy_bis[4]")
               .Define("jetconstituents_antikt4_energy_6_bis", "jetconstituents_antikt4_energy_bis[5]")

               
               .Define("the2Z_antikt_bis", "ReconstructedParticle::resoantikt(jets_antikt_e4_bis, jets_antikt_px4_bis, jets_antikt_py4_bis ,jets_antikt_pz4_bis, N_jets_antikt4_Ab5_bis, jets_antikt_theta4_bis)")

               .Define("firstZ_antikt_bis", "the2Z_antikt_bis.Z1")
               .Define("firstZ_m_antikt_bis","firstZ_antikt_bis.M()")

               .Define("secondZ_antikt_bis", "the2Z_antikt_bis.Z2")
               .Define("secondZ_m_antikt_bis","secondZ_antikt_bis.M()")

               .Define("firstZ_antikt_1_2_bis", "the2Z_antikt_bis.Z1_1_2")
               .Define("firstZ_m_antikt_1_2_bis", "firstZ_antikt_1_2_bis.M()")

               .Define("firstZ_antikt_reco_bis", "the2Z_antikt_bis.Z1_reco")
               .Define("firstZ_m_antikt_reco_bis", "firstZ_antikt_reco_bis.M()")
               

            

                #fin

                #Durham anti-kt R=0.8

               .Define("N_jets_antikt4_Ab5_ter", "ReconstructedParticle::countNjets(jets_antikt_e4_ter, 5)")
               .Define("N_jets_antikt4_Ab10_ter", "ReconstructedParticle::countNjets(jets_antikt_e4_ter, 10)")

               .Define("the2Z_antikt_ter", "ReconstructedParticle::resoantikt(jets_antikt_e4_ter, jets_antikt_px4_ter, jets_antikt_py4_ter ,jets_antikt_pz4_ter, N_jets_antikt4_Ab5_ter, jets_antikt_theta4_ter)")

               .Define("firstZ_antikt_ter", "the2Z_antikt_ter.Z1")
               .Define("firstZ_m_antikt_ter","firstZ_antikt_ter.M()")

               .Define("secondZ_antikt_ter", "the2Z_antikt_ter.Z2")
               .Define("secondZ_m_antikt_ter","secondZ_antikt_ter.M()")

               .Define("firstZ_antikt_1_2_ter", "the2Z_antikt_ter.Z1_1_2")
               .Define("firstZ_m_antikt_1_2_ter", "firstZ_antikt_1_2_ter.M()")

               .Define("firstZ_antikt_reco_ter", "the2Z_antikt_ter.Z1_reco")
               .Define("firstZ_m_antikt_reco_ter", "firstZ_antikt_reco_ter.M()")
               
                #fin

                #retour Durham kt N=4

               .Define("jetconstituents_4_theta_1", "jetconstituents_4_theta[0]")
               .Define("jetconstituents_4_theta_2", "jetconstituents_4_theta[1]")
               .Define("jetconstituents_4_theta_3", "jetconstituents_4_theta[2]")
               .Define("jetconstituents_4_theta_4", "jetconstituents_4_theta[3]")

               .Define("jetconstituents_4_phi_1", "jetconstituents_4_phi[0]")
               .Define("jetconstituents_4_phi_2", "jetconstituents_4_phi[1]")
               .Define("jetconstituents_4_phi_3", "jetconstituents_4_phi[2]")
               .Define("jetconstituents_4_phi_4", "jetconstituents_4_phi[3]")

               .Define("jetconstituents_4_energy_1", "jetconstituents_4_energy[0]")
               .Define("jetconstituents_4_energy_2", "jetconstituents_4_energy[1]")
               .Define("jetconstituents_4_energy_3", "jetconstituents_4_energy[2]")
               .Define("jetconstituents_4_energy_4", "jetconstituents_4_energy[3]")

               .Define("firstZ_jets_eta_vector", "the2Z.etaZ1")
               .Define("deltaetaZ1", "abs(firstZ_jets_eta_vector[0] - firstZ_jets_eta_vector[1])")
               .Define("secondZ_jets_eta_vector", "the2Z.etaZ2")
               .Define("deltaetaZ2", "abs(secondZ_jets_eta_vector[0] - secondZ_jets_eta_vector[1])")

               .Define("firstZ", "the2Z.Z1")
               .Define("firstZ_m", "firstZ.M()")
               .Define("firstZ_px","firstZ.Px()")
               .Define("firstZ_py","firstZ.Py()")
               .Define("firstZ_pz","firstZ.Pz()")
               .Define("firstZ_p", "firstZ.P()")
               .Define("firstZ_pt","firstZ.Pt()")
               .Define("firstZ_theta", "firstZ.Theta()")
               .Define("firstZ_eta", "firstZ.Eta()")

               .Define("indicesofthejets", "the2Z.jetmember")
               .Define("firstZ_firstjet", "indicesofthejets[0]")
               .Define("firstZ_secondjet", "indicesofthejets[1]")
               .Define("secondZ_firstjet", "indicesofthejets[2]")
               .Define("secondZ_secondjet", "indicesofthejets[3]")


#flavour of the jets of the 2 Z have been stored in the [2] for Z1 and [3] for Z2 in the output of myresobuilder 
#they are tlorentzvectors (bizarre), in the Px entry is is the flavor of the 1st jet and Py the other jet

               .Define("firstZ_flavour_vector", "the2Z.flav1")
               .Define("firstZ_part1_flavour", "firstZ_flavour_vector[0]")
               .Define("firstZ_part2_flavour", "firstZ_flavour_vector[1]")
              
               .Define("firstZ_part1_flavour_gm", "firstZ_flavour_vector[2]")
               .Define("firstZ_part2_flavour_gm", "firstZ_flavour_vector[3]")
            
               
               .Define("secondZ", "the2Z.Z2")
               .Define("secondZ_m", "secondZ.M()")
               .Define("secondZ_px","secondZ.Px()")
               .Define("secondZ_py","secondZ.Py()")
               .Define("secondZ_pz","secondZ.Pz()")
               .Define("secondZ_p", "secondZ.P()")
               .Define("secondZ_pt","secondZ.Pt()")
               .Define("secondZ_theta", "secondZ.Theta()")
               .Define("secondZ_eta", "secondZ.Eta()")

               .Define("secondZ_flavour_vector", "the2Z.flav2") 
               .Define("secondZ_part1_flavour", "secondZ_flavour_vector[0]")
               .Define("secondZ_part2_flavour", "secondZ_flavour_vector[1]")

 	       .Define("secondZ_part1_flavour_gm", "secondZ_flavour_vector[2]")
               .Define("secondZ_part2_flavour_gm", "secondZ_flavour_vector[3]")



               #angular differences between the paired jets

               .Define("angulardiff", "the2Z.angulardiff")
               .Define("Dtheta_Z1", "angulardiff[0]")
               .Define("Dphi_Z1", "angulardiff[1]") 
               .Define("Dtheta_Z2", "angulardiff[2]")
               .Define("Dphi_Z2", "angulardiff[3]") 

               #min and second min angular difference between the jets 
               .Define("mintheta", "angulardiff[4]") 
               .Define("secondmintheta", "angulardiff[5]") 


#TLorentzVector du higgs que l'on reconstruit en sommant les deux Z de son état final 
               
               .Define("higgs_4", "firstZ + secondZ")
               .Define("higgs_4_m", "higgs_4.M()")
               .Define("higgs_4_px", "higgs_4.Px()")
               .Define("higgs_4_py", "higgs_4.Py()")
               .Define("higgs_4_pz", "higgs_4.Pz()")
               .Define("higgs_4_p", "higgs_4.P()")
               .Define("deltaetaZZ", "abs(firstZ_eta - secondZ_eta)")

               #.Define("diffmasshiggs", "abs( Zcand_recoil_m - higgs_4_m )") 

#delta theta entre chaque Z, Z étant le leptonique, Z1 le ON SHELL du Higgs et Z2 le OFF SHELL du Higgs

               .Define("Dtheta_ZZ1", "abs(Zcand_theta - firstZ_theta)")
               .Define("Dtheta_ZZ2", "abs(Zcand_theta - secondZ_theta)") 
               .Define("Dtheta_Z1Z2", "abs( firstZ_theta - secondZ_theta)")

               .Define("flavourscore", "ReconstructedParticle::sameflavour()(firstZ_part1_flavour, firstZ_part2_flavour, secondZ_part1_flavour, secondZ_part2_flavour)") 
               .Define("flavourscoreZ1", "flavourscore[0]")
               .Define("flavourscoreZ2", "flavourscore[1]")

               .Define("flavourscoregm", "ReconstructedParticle::sameflavour()(firstZ_part1_flavour_gm, firstZ_part2_flavour_gm, secondZ_part1_flavour_gm, secondZ_part2_flavour_gm)") 
               .Define("flavourscoreZ1_gm", "flavourscoregm[0]")
               .Define("flavourscoreZ2_gm", "flavourscoregm[1]")

               
               #durham trie naturellement les jets dans l'ordre décroissant de p

               .Define("jet1_p", "jets_p4[0]")
               .Define("jet2_p", "jets_p4[1]")
               .Define("jet3_p", "jets_p4[2]")
               .Define("jet4_p", "jets_p4[3]")

               .Define("jet1_pt", "jets_pt4[0]")
               .Define("jet2_pt", "jets_pt4[1]")
               .Define("jet3_pt", "jets_pt4[2]")
               .Define("jet4_pt", "jets_pt4[3]")

               .Define("jet1_e", "jets_e4[0]")
               .Define("jet2_e", "jets_e4[1]")
               .Define("jet3_e", "jets_e4[2]")
               .Define("jet4_e", "jets_e4[3]")

               .Define("jet1_px", "jets_px4[0]")
               .Define("jet2_px", "jets_px4[1]")
               .Define("jet3_px", "jets_px4[2]")
               .Define("jet4_px", "jets_px4[3]")
               
               .Define("jet1_py", "jets_py4[0]")
               .Define("jet2_py", "jets_py4[1]")
               .Define("jet3_py", "jets_py4[2]")
               .Define("jet4_py", "jets_py4[3]")

               .Define("jet1_pz", "jets_pz4[0]")
               .Define("jet2_pz", "jets_pz4[1]")
               .Define("jet3_pz", "jets_pz4[2]")
               .Define("jet4_pz", "jets_pz4[3]")

               .Define("jet1_theta", "jets_theta4[0]")
               .Define("jet2_theta", "jets_theta4[1]")
               .Define("jet3_theta", "jets_theta4[2]")
               .Define("jet4_theta", "jets_theta4[3]")
               
               .Define("jet1_phi", "jets_phi4[0]")
               .Define("jet2_phi", "jets_phi4[1]")
               .Define("jet3_phi", "jets_phi4[2]")
               .Define("jet4_phi", "jets_phi4[3]")


               

               #2jets de Durham

               #on met les 2 jets dans le bon ordre de masse, cad le plus proche de la masse du Z en premier pour avoir Z1 et Z2 en second

               .Define("jets_2_ordered", "ReconstructedParticle::mass_order(91)(jets_m2, jets_px2, jets_py2, jets_pz2, jets_e2)")
               .Define("jet1_2", "jets_2_ordered[0]")
               .Define("jet1_2_m", "jet1_2.M()")
               .Define("jet1_2_px","jet1_2.Px()")
               .Define("jet1_2_py","jet1_2.Py()")
               .Define("jet1_2_pz","jet1_2.Pz()")
               .Define("jet1_2_e","jet1_2.E()")

               .Define("jet2_2", "jets_2_ordered[1]")
               .Define("jet2_2_m", "jet2_2.M()")
               .Define("jet2_2_px","jet2_2.Px()")
               .Define("jet2_2_py","jet2_2.Py()")
               .Define("jet2_2_pz","jet2_2.Pz()")
               .Define("jet2_2_e","jet2_2.E()")

               .Define("higgs_2", "jet1_2 + jet2_2")
               .Define("higgs_2_m", "higgs_2.M()")            

               #test efficacité de durham 4 pour reconstruire correctement les jets 

               #.Define("diffZ14", "abs(inv_mass_Z1 - firstZ_m)")
               #.Define("diffZ24", "abs(inv_mass_Z2 - secondZ_m)")

               
               #test efficacité de durham 2 pour reconstruire les bons Z

               #.Define("diffZ12", "abs( inv_mass_Z1 - jet1_2_m )")
               #.Define("diffZ22", "abs( inv_mass_Z2 - jet2_2_m )")

               .Define("thetaphi", "MCParticle::fill_thetaphi_ZHdecay(ZH_decay)")
               .Define("Hpart_theta", "thetaphi.H_theta")
               .Define("Hpart_phi", "thetaphi.H_phi")
               .Define("Hpart_energy", "thetaphi.H_energy")
               .Define("Zpart_theta", "thetaphi.Z_theta")
               .Define("Zpart_phi", "thetaphi.Z_phi")
               .Define("Zpart_energy", "thetaphi.Z_energy")
               
               .Define("truth_Zq1_theta", "thetaphi.truth_Zq1_theta")
               .Define("truth_Zq1_phi", "thetaphi.truth_Zq1_phi")
               
               .Define("truth_Zq2_theta", "thetaphi.truth_Zq2_theta")
               .Define("truth_Zq2_phi", "thetaphi.truth_Zq2_phi")
               
               .Define("truth_Hq1_theta", "thetaphi.truth_Hq1_theta")
               .Define("truth_Hq1_phi", "thetaphi.truth_Hq1_phi")
               
               .Define("truth_Hq2_theta", "thetaphi.truth_Hq2_theta")
               .Define("truth_Hq2_phi", "thetaphi.truth_Hq2_phi")


        
               )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list.
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

            "N_zed_leptonic",
            "N_selected_leptons",
            
            "Zcand_m",
            "Zcand_recoil_m",
            "Zcand_pt",
            "Zcand_px",
            "Zcand_py",
            "Zcand_pz",
            "Zcand_p",
            "Zcand_q",
            "Zcand_phi",
            "Zcand_theta",
            "Zcand_cos",
            "Zcand_y",
            "Zcand_eta",

            "firstZ_m",
            "firstZ_px",
            "firstZ_py",
            "firstZ_pz",
            "firstZ_p",
            "firstZ_pt",
            "firstZ_theta",
            "firstZ_eta",

            "firstZ_part1_flavour",
            "firstZ_part2_flavour",
            "firstZ_part1_flavour_gm",
            "firstZ_part2_flavour_gm",

            
            "secondZ_m",
            "secondZ_px",
            "secondZ_py",
            "secondZ_pz",
            "secondZ_p",
            "secondZ_pt",
            "secondZ_theta",  
            "secondZ_eta",
            
            "secondZ_part1_flavour",
            "secondZ_part2_flavour",
            "secondZ_part1_flavour_gm",
            "secondZ_part2_flavour_gm",

            "jet1_p",
            "jet1_pt",
            "jet1_e",
            "jet1_px",
            "jet1_py",
            "jet1_pz",

            "jet2_p",
            "jet2_pt",
            "jet2_e",
            "jet2_px",
            "jet2_py",
            "jet2_pz",

            "jet3_p",
            "jet3_pt",
            "jet3_e",
            "jet3_px",
            "jet3_py",
            "jet3_pz",

            "jet4_p",
            "jet4_pt",
            "jet4_e",
            "jet4_px",
            "jet4_py",
            "jet4_pz",

            "higgs_4_m",
            "higgs_4_p",
            "higgs_4_px",
            "higgs_4_py",
            "higgs_4_pz",
            "deltaetaZZ",
            "deltaetaZ1",
            "deltaetaZ2",

            "jet1_2_m",
            "jet2_2_m",
            "higgs_2_m",

            #"diffZ14",
            #"diffZ24",
            
            #"diffZ12",
            #"diffZ22",
            
            "emiss",
            "pzmiss",
            "etmiss",
            
            
            "dmerge_2_45",
            "dmerge_2_34",
            "dmerge_2_23",
            "dmerge_2_12",

            "dmerge_4_45",
            "dmerge_4_34",
            "dmerge_4_23",
            "dmerge_4_12",

        
            "Dtheta_ZZ1",
            "Dtheta_ZZ2", 
            "Dtheta_Z1Z2",

            "Dtheta_Z1",
            "Dphi_Z1",
            "Dtheta_Z2",
            "Dphi_Z2",
            "mintheta",
            "secondmintheta", 

            "flavourscoreZ1",
            "flavourscoreZ2",
            "flavourscoreZ1_gm",
            "flavourscoreZ2_gm",

            "Hpart_theta",
            "Hpart_phi",
            "Zpart_theta",
            "Zpart_phi",
            
            "truth_Zq1_theta",
            "truth_Zq1_phi",
            "truth_Zq2_theta",
            "truth_Zq2_phi",
            
            "truth_Hq1_theta",
            "truth_Hq1_phi",
            "truth_Hq2_theta",
            "truth_Hq2_phi",

            "jetconstituents_4_1",
            "jetconstituents_4_2",
            "jetconstituents_4_3",
            "jetconstituents_4_4",

            "jet1_theta",
            "jet2_theta",
            "jet3_theta",
            "jet4_theta",

            "jet1_phi",
            "jet2_phi",
            "jet3_phi",
            "jet4_phi",

            "jetconstituents_4_theta",
            "jetconstituents_4_phi",
            "jetconstituents_4_energy",

            "jetconstituents_4_theta_1",
            "jetconstituents_4_theta_2",
            "jetconstituents_4_theta_3",
            "jetconstituents_4_theta_4",

            "jetconstituents_4_phi_1",
            "jetconstituents_4_phi_2",
            "jetconstituents_4_phi_3",
            "jetconstituents_4_phi_4",

            "jetconstituents_4_energy_1",
            "jetconstituents_4_energy_2",
            "jetconstituents_4_energy_3",
            "jetconstituents_4_energy_4",

            
            "Hpart_energy",
            "Zpart_energy",

            "firstZ_firstjet",
            "firstZ_secondjet",
            "secondZ_firstjet",
            "secondZ_secondjet",

            #durham antikt R=0.4

            "N_jets_antikt4",
            "N_jets_antikt4_Ab5",
            "N_jets_antikt4_Ab10",

            "jets_antikt_e4",
            "jet1_antikt4_energy",
            "jet2_antikt4_energy",
            "jet3_antikt4_energy",
            "jet4_antikt4_energy",
            "jet5_antikt4_energy",
            "jet6_antikt4_energy",
            "jet7_antikt4_energy",
            "jet8_antikt4_energy",
            "jet9_antikt4_energy",
            "jet10_antikt4_energy",
            "jet11_antikt4_energy",
            "jet12_antikt4_energy",

            "jetconstituents_antikt4",
            "jetconstituents_antikt4_1",
            "jetconstituents_antikt4_2",
            "jetconstituents_antikt4_3",
            "jetconstituents_antikt4_4",

            "jetconstituents_antikt4_theta",
            "jetconstituents_antikt4_theta_1",
            "jetconstituents_antikt4_theta_2",
            "jetconstituents_antikt4_theta_3",
            "jetconstituents_antikt4_theta_4",

            "jetconstituents_antikt4_phi",
            "jetconstituents_antikt4_phi_1",
            "jetconstituents_antikt4_phi_2",
            "jetconstituents_antikt4_phi_3",
            "jetconstituents_antikt4_phi_4",

            "jetconstituents_antikt4_energy",
            "jetconstituents_antikt4_energy_1",
            "jetconstituents_antikt4_energy_2",
            "jetconstituents_antikt4_energy_3",
            "jetconstituents_antikt4_energy_4",

            "dmerge_antikt4_45",
            "dmerge_antikt4_34",
            "dmerge_antikt4_23",
            "dmerge_antikt4_12",
        

            "firstZ_m_antikt",
            "secondZ_m_antikt",
            "firstZ_m_antikt_1_2",
            "firstZ_m_antikt_reco",

            #durham antikt R=0.6

            "N_jets_antikt4_bis",
            "N_jets_antikt4_Ab5_bis",
            "N_jets_antikt4_Ab10_bis",

            "jets_antikt_e4_bis",
            "jet1_antikt4_energy_bis",
            "jet2_antikt4_energy_bis",
            "jet3_antikt4_energy_bis",
            "jet4_antikt4_energy_bis",
            "jet5_antikt4_energy_bis",
            "jet6_antikt4_energy_bis",
            "jet7_antikt4_energy_bis",
            "jet8_antikt4_energy_bis",
            "jet9_antikt4_energy_bis",
            "jet10_antikt4_energy_bis",
            "jet11_antikt4_energy_bis",
            "jet12_antikt4_energy_bis",

            "jetconstituents_antikt4_bis",
            "jetconstituents_antikt4_1_bis",
            "jetconstituents_antikt4_2_bis",
            "jetconstituents_antikt4_3_bis",
            "jetconstituents_antikt4_4_bis",
            "jetconstituents_antikt4_5_bis",
            "jetconstituents_antikt4_6_bis",

            "jetconstituents_antikt4_theta_bis",
            "jetconstituents_antikt4_theta_1_bis",
            "jetconstituents_antikt4_theta_2_bis",
            "jetconstituents_antikt4_theta_3_bis",
            "jetconstituents_antikt4_theta_4_bis",
            "jetconstituents_antikt4_theta_5_bis",
            "jetconstituents_antikt4_theta_6_bis",

            "jetconstituents_antikt4_phi_bis",
            "jetconstituents_antikt4_phi_1_bis",
            "jetconstituents_antikt4_phi_2_bis",
            "jetconstituents_antikt4_phi_3_bis",
            "jetconstituents_antikt4_phi_4_bis",
            "jetconstituents_antikt4_phi_5_bis",
            "jetconstituents_antikt4_phi_6_bis",

            "jetconstituents_antikt4_energy_bis",
            "jetconstituents_antikt4_energy_1_bis",
            "jetconstituents_antikt4_energy_2_bis",
            "jetconstituents_antikt4_energy_3_bis",
            "jetconstituents_antikt4_energy_4_bis",
            "jetconstituents_antikt4_energy_5_bis",
            "jetconstituents_antikt4_energy_6_bis",

            "dmerge_antikt4_45_bis",
            "dmerge_antikt4_34_bis",
            "dmerge_antikt4_23_bis",
            "dmerge_antikt4_12_bis",

            "firstZ_m_antikt_bis",
            "secondZ_m_antikt_bis",
            "firstZ_m_antikt_1_2_bis",
            "firstZ_m_antikt_reco_bis",
            
            "N_jets_antikt4_ter",
            "N_jets_antikt4_Ab5_ter",
            "N_jets_antikt4_Ab10_ter",

            "firstZ_m_antikt_ter",
            "secondZ_m_antikt_ter",
            "firstZ_m_antikt_1_2_ter",
            "firstZ_m_antikt_reco_ter",

            
            "truth_Z_theta",
            "truth_Z_phi",
            
            "truth_H_theta",
            "truth_H_phi",


            "N_jets_4" #test pour durham 4 

            

            
        ]
        return branchList
