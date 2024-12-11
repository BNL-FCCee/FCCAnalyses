inputDir    = "./2DPlotting/stage1"

#Optional: output directory, default is local dir
outputDir    = "./2DPlotting/stage2"
# outputDir   = "outputs/fccee/higgs/mH-recoil/hzz/stage2/4jets/Zleptonic/"

processList = {}

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
#maintenant, on est plus que sur des évènements où le Z seul est en leptons et le higgs en deux Z eux meme en 2 jets (4jets)
#Now we are on events where only the Z is of leptons and the Higgs of 2 Z which each have 2 jets which gives 4 jets
               .Filter("size(jets_e4)>0") 
               #.Define Z candidate mass
               .Define("Zcand_m","zed_leptonic_m[0]")
               #Define Z candidate recoil mass
               .Define("Zcand_recoil_m","zed_leptonic_recoil_m[0]")
               #Define Z candidate pt
               .Define("Zcand_pt","zed_leptonic_pt[0]")
                #px
               .Define("Zcand_px","zed_leptonic_px[0]")
                #py
               .Define("Zcand_py","zed_leptonic_py[0]")
               #pz
               .Define("Zcand_pz","zed_leptonic_pz[0]")
               #p
               .Define("Zcand_p","zed_leptonic_p[0]")
               #charge
               .Define("Zcand_q","zed_leptonic_charge[0]")
               #phi
               .Define("Zcand_phi", "zed_leptonic_phi[0]")
               #theta
               .Define("Zcand_theta", "zed_leptonic_theta[0]")
               #cos
               .Define("Zcand_cos","zed_leptonic_cos[0]")
               #y
               .Define("Zcand_y", "zed_leptonic_y[0]")
               #eta
               .Define("Zcand_eta", "zed_leptonic_eta[0]")

                #resonance Builder: we regroup the jets by 2

                .Define("the2Z", "ReconstructedParticle::myresoBuilder(jets_e4, jets_px4, jets_py4, jets_pz4, jets_ee_genkt_flavour4, jets_ee_flavour4, jets_eta4, jets_theta4, jets_phi4)")

                 #antikt R=0.4

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

               #Reconstructing of Zs from the jets
               

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


                # Durham kt N=4

               .Define("jetconstituents_4_1", "jetconstituents_4[0]")
               .Define("jetconstituents_4_2", "jetconstituents_4[1]")
               .Define("jetconstituents_4_3", "jetconstituents_4[2]")
               .Define("jetconstituents_4_4", "jetconstituents_4[3]")

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

               
             #durham naturally sorts the jets in decreasing order of p
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

            #antikt R=0.4

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
            

            "truth_Z_theta",
            "truth_Z_phi",
            "truth_H_theta",
            "truth_H_phi",

            "N_jets_4" #test pour durham 4 


            
        ]
        return branchList
