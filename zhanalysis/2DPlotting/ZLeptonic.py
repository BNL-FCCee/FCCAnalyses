inputDir    = "./2DPlotting/antikt/stage1"

#Optional: output directory, default is local dir
outputDir    = "./2DPlotting/antikt/stage2"

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
               
        #Durham Kt       
               #defining columns for the jet constituent arrays for each jet
               
               #.Define("jetconstituents_kt4_1", "jetconstituents_kt4[0]")
               #.Define("jetconstituents_kt4_2", "jetconstituents_kt4[1]")
               #.Define("jetconstituents_kt4_3", "jetconstituents_kt4[2]")
               #.Define("jetconstituents_kt4_4", "jetconstituents_kt4[3]")

               #jet consituents theta 
               .Define("jetconstituents_kt4_theta_1", "jetconstituents_kt4_theta[0]")
               .Define("jetconstituents_kt4_theta_2", "jetconstituents_kt4_theta[1]")
               .Define("jetconstituents_kt4_theta_3", "jetconstituents_kt4_theta[2]")
               .Define("jetconstituents_kt4_theta_4", "jetconstituents_kt4_theta[3]")

                #jet consituents phi
               .Define("jetconstituents_kt4_phi_1", "jetconstituents_kt4_phi[0]")
               .Define("jetconstituents_kt4_phi_2", "jetconstituents_kt4_phi[1]")
               .Define("jetconstituents_kt4_phi_3", "jetconstituents_kt4_phi[2]")
               .Define("jetconstituents_kt4_phi_4", "jetconstituents_kt4_phi[3]")

                #jet consituents phi
               .Define("jetconstituents_kt4_energy_1", "jetconstituents_kt4_energy[0]")
               .Define("jetconstituents_kt4_energy_2", "jetconstituents_kt4_energy[1]")
               .Define("jetconstituents_kt4_energy_3", "jetconstituents_kt4_energy[2]")
               .Define("jetconstituents_kt4_energy_4", "jetconstituents_kt4_energy[3]")


            #construct the 2 Z
               .Define("the2Z", "ReconstructedParticle::myresoBuilder(jets_kt_e4, jets_kt_px4, jets_kt_py4 ,jets_kt_pz4, jets_kt_flavour_4, jets_kt_flavour_4, jets_kt_eta4, jets_kt_theta4, jets_kt_phi4)")
              
              #first Z and mass
               .Define("firstZ_kt", "the2Z.Z1")
               .Define("firstZ_kt_m", "firstZ_kt.M()")
              #second Z and mass
               .Define("secondZ_kt", "the2Z.Z2")
               .Define("secondZ_kt_m", "secondZ_kt.M()")
            
               .Define("higgs_kt4", "firstZ_kt + secondZ_kt")
               .Define("higgs_kt4_m", "higgs_kt4.M()")

               .Define("indicesofthejets", "the2Z.jetmember")
               .Define("firstZ_firstjet", "indicesofthejets[0]")
               .Define("firstZ_secondjet", "indicesofthejets[1]")
               .Define("secondZ_firstjet", "indicesofthejets[2]")
               .Define("secondZ_secondjet", "indicesofthejets[3]")                  



#flavour of the jets of the 2 Z have been stored in the [2] for Z1 and [3] for Z2 in the output of myresobuilder 
#they are tlorentzvectors (bizarre), in the Px entry is is the flavor of the 1st jet and Py the other jet

            
        #Antikt


               .Define("N_jets_Ab5", "ReconstructedParticle::countNjets(jets_antikt_e4, 5)")

               .Define("the2Z_antikt", "ReconstructedParticle::resoantikt(jets_antikt_e4, jets_antikt_px4, jets_antikt_py4 ,jets_antikt_pz4, N_jets_Ab5, jets_antikt_theta4)")

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


               #jet consituents theta 
               .Define("jetconstituents_antikt4_theta_1", "jetconstituents_antikt4_theta[0]")
               .Define("jetconstituents_antikt4_theta_2", "jetconstituents_antikt4_theta[1]")
               .Define("jetconstituents_antikt4_theta_3", "jetconstituents_antikt4_theta[2]")
               .Define("jetconstituents_antikt4_theta_4", "jetconstituents_antikt4_theta[3]")

                #jet consituents phi
               .Define("jetconstituents_antikt4_phi_1", "jetconstituents_antikt4_phi[0]")
               .Define("jetconstituents_antikt4_phi_2", "jetconstituents_antikt4_phi[1]")
               .Define("jetconstituents_antikt4_phi_3", "jetconstituents_antikt4_phi[2]")
               .Define("jetconstituents_antikt4_phi_4", "jetconstituents_antikt4_phi[3]")

                #jet consituents phi
               .Define("jetconstituents_antikt4_energy_1", "jetconstituents_antikt4_energy[0]")
               .Define("jetconstituents_antikt4_energy_2", "jetconstituents_antikt4_energy[1]")
               .Define("jetconstituents_antikt4_energy_3", "jetconstituents_antikt4_energy[2]")
               .Define("jetconstituents_antikt4_energy_4", "jetconstituents_antikt4_energy[3]")


               #methode : Z1 : 1+2+5+6+... ; Z2 = 3 + 4

              #first Z and mass
               .Define("firstZ_antikt", "the2Z_antikt.Z1")
               .Define("firstZ_antikt_m", "firstZ_antikt.M()")
              #second Z and mass
               .Define("secondZ_antikt", "the2Z.Z2")
               .Define("secondZ_antikt_m", "secondZ_antikt.M()")
            
               .Define("higgs_antikt4", "firstZ_antikt + secondZ_antikt")
               .Define("higgs_antikt4_m", "higgs_antikt4.M()")

        
               #methode : Z1 = 1+2

               .Define("firstZ_antikt_1_2", "the2Z_antikt.Z1_1_2")
               .Define("firstZ_antikt_m_1_2", "firstZ_antikt_1_2.M()")

               #methode : Z1 = 1+2 + tout ce qui est + proche de 1 ou 2 que de 3 ou 4 ; Z2 = 3+4+ ce qui n'est pas avec Z1

               .Define("firstZ_antikt_reco", "the2Z_antikt.Z1_reco")
               .Define("firstZ_antikt_m_reco", "firstZ_antikt_reco.M()")


            #truth particles

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
            "Hpart_phi",
            "Hpart_theta",
            "Hpart_energy",
            "Zpart_phi",
            "Zpart_theta",
            "Zpart_energy",
            
            "jetconstituents_antikt4_theta_1",
            "jetconstituents_antikt4_theta_2",
            "jetconstituents_antikt4_theta_3",
            "jetconstituents_antikt4_theta_4",

            "jetconstituents_antikt4_energy_1",
            "jetconstituents_antikt4_energy_2",
            "jetconstituents_antikt4_energy_3",
            "jetconstituents_antikt4_energy_4",

            "jetconstituents_antikt4_phi_1",
            "jetconstituents_antikt4_phi_2",
            "jetconstituents_antikt4_phi_3",
            "jetconstituents_antikt4_phi_4",

            "firstZ_antikt_m",
            "secondZ_antikt_m",
            "higgs_antikt4_m",
            
            "firstZ_antikt_1_2",
            "firstZ_antikt_m_1_2",
            "firstZ_antikt_reco",
            "firstZ_antikt_m_reco",
    
            "jetconstituents_kt4_theta_1",
            "jetconstituents_kt4_theta_2",
            "jetconstituents_kt4_theta_3",
            "jetconstituents_kt4_theta_4",

            "jetconstituents_kt4_energy_1",
            "jetconstituents_kt4_energy_2",
            "jetconstituents_kt4_energy_3",
            "jetconstituents_kt4_energy_4",

            "jetconstituents_kt4_phi_1",
            "jetconstituents_kt4_phi_2",
            "jetconstituents_kt4_phi_3",
            "jetconstituents_kt4_phi_4",

            "firstZ_kt_m",
            "secondZ_kt_m",
            "higgs_kt4_m",


            "truth_Zq1_theta",
            "truth_Zq1_phi",
            "truth_Zq2_theta",
            "truth_Zq2_phi",

            "truth_H_theta",
            "truth_H_phi",
            "truth_Hq1_theta",
            "truth_Hq1_phi",
            "truth_Hq2_theta",
            "truth_Hq2_phi"
    
           # ,"firstZ_firstjet",
           # "firstZ_secondjet",
           # "secondZ_firstjet",
           # "secondZ_secondjet",

        
        ]
        
        return branchList
