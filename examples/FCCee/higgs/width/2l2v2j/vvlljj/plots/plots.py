import ROOT

# global parameters
intLumi        =  5.0e+06 #in pb-1
intLumiLabel   = "L = 5 ab^{-1}"
ana_tex        = "ZH #rightarrow ZZZ #rightarrow #nu#nu + #mu^{+}#mu^{-}/e^{+}e^{-} jj"
delphesVersion = "3.4.2"
energy         =  240.0
collider       = "FCC-ee"
inputDir       = "outputs/fccee/higgs/mH-recoil/neutrinos/vvlljj/final/"
formats        = ['png','pdf']
yaxis          = ['lin','log']
stacksig       = ['stack','nostack']
outdir         = "outputs/fccee/higgs/mH-recoil/neutrinos/vvlljj/plots/"

variables = ['mz','mz_zoom', 'theta', 'phi', 'cos' ,'eta', 'y','leptonic_recoil_m','leptonic_recoil_m_zoom','leptonic_recoil_m_zoom2','leptonic_recoil_m_zoom6', 'selected_muons_e', 'selected_electrons_e', 'selected_leptons_e', 'selected_muons_p', 'selected_electrons_p' , 'selected_leptons_p', 'selected_muons_px', 'selected_electrons_px' , 'selected_leptons_px', 'selected_muons_py', 'selected_electrons_py' , 'selected_leptons_py', 'selected_muons_pz', 'selected_electrons_pz' , 'selected_leptons_pz', 'selected_muons_pt', 'selected_electrons_pt' , 'selected_leptons_pt', 'N_selected_leptons', 'N_zed_leptonic', 'jet1_p', 'jet2_p', 'jet1_px', 'jet2_px', 'jet1_py', 'jet2_py' , 'jet1_pz', 'jet2_pz', 'jet1_e', 'jet2_e', 'jet1_pt', 'jet2_pt', 'jet1_theta', 'jet2_theta', 'jet1_Nconst', 'jet2_Nconst', 'diffthetajets', 'secondZ_m', 'secondZ_m_zoom1', 'secondZ_m_zoom2', 'secondZ_p', 'secondZ_px', 'secondZ_py', 'secondZ_pz', 'secondZ_pt', 'dmerge_2_12', 'dmerge_2_23', 'dmerge_2_34', 'dmerge_2_45', 'emiss', 'pxmiss', 'pymiss', 'pzmiss', 'etmiss', 'etmiss_zoom', 'visible_mass_predef', 'visible_mass_predef_zoom', 'N_LooseLeptons', 'LooseLeptons_pt', 'LooseLeptons_p', 'LooseLeptons_phi', 'LooseLeptons_theta', 'minNconst', 'maxNconst', 'meanNconst', 'minNconst_3', 'maxNconst_3', 'meanNconst_3', 'missing_theta', 'angle_miss_jet1', 'angle_miss_jet2', 'min_angle_miss_jet', 'max_angle_miss_jet', 'N_LooseLeptons_2','dmerge_2_12_zoom', 'visible_mass_predef_centered']

###Dictonnary with the analysis name as a key, and the list of selections to be plotted for this analysis. The name of the selections should be the same than in the final selection
selections = {}
selections['vvlljj']   = ["sel0","sel1","sel2","sel3","sel4","sel5","sel6", "sel7","sel8"]
#,"sel2","sel3", "sel4", "sel5", "sel6", "sel7","sel8","sel9","sel10","sel11"]


extralabel = {}
extralabel['sel0'] = "no cut"
extralabel['sel1'] = "Nb lep"
extralabel['sel2'] = "Nb lep + emiss"
extralabel['sel3'] = "Nb lep + emiss + VisMass"
extralabel['sel4'] = "Nb lep + emiss + VisMass + mjj"
extralabel['sel5'] = "Nb lep + emiss + VisMass + mjj + MeanNconst"
extralabel['sel6'] = "Nb lep + emiss + VisMass + mjj + MeanNconst + mll"
extralabel['sel7'] = "Nb lep + emiss + VisMass + mjj + MeanNconst + mll + etmiss"
extralabel['sel8'] = "Final selection"
extralabel['sel9'] = "Nb lep + mll + recoil + MeanNconst + mjj + cosMissAng + minAngle + NLeptons(p>2) + emiss"
extralabel['sel10'] = "Region 1, dmerge_23 < 60"
extralabel['sel11'] = "Region 2, dmerge_23 > 60"





colors = {}
colors['Signal,HZZ'] = ROOT.kBlue
colors['llH,HZZ'] = ROOT.kRed
colors['WW'] = ROOT.kOrange+6
colors['ZZ'] = ROOT.kOrange-2
colors['Htautau'] = ROOT.kCyan+2
colors['HWW'] = ROOT.kGreen
colors['Hbb'] = ROOT.kPink-9
colors['H->other decays'] = ROOT.kViolet

#colors['VV'] = ROOT.kGreen+3

plots = {}
plots['vvlljj'] = {'signal':{'Signal,HZZ':['wzp6_ee_nunuH_HZZ_ecm240']},
                 'backgrounds':{'llH,HZZ':['wzp6_ee_mumuH_HZZ_ecm240','wzp6_ee_eeH_HZZ_ecm240'],

                                   'WW':['p8_ee_WW_ecm240'],
                                   'ZZ':['p8_ee_ZZ_ecm240'], 

                                   'Htautau':['wzp6_ee_mumuH_Htautau_ecm240', 'wzp6_ee_nunuH_Htautau_ecm240', 'wzp6_ee_eeH_Htautau_ecm240'],
                                   
                                   'HWW':['wzp6_ee_mumuH_HWW_ecm240', 'wzp6_ee_nunuH_HWW_ecm240', 'wzp6_ee_eeH_HWW_ecm240'],
                                   
                                   'Hbb':['wzp6_ee_mumuH_Hbb_ecm240', 'wzp6_ee_nunuH_Hbb_ecm240', 'wzp6_ee_eeH_Hbb_ecm240'],
                                   
                                   'H->other decays':['wzp6_ee_mumuH_Hcc_ecm240', 'wzp6_ee_nunuH_Hcc_ecm240', 'wzp6_ee_eeH_Hcc_ecm240', 'wzp6_ee_mumuH_Hgg_ecm240', 'wzp6_ee_nunuH_Hgg_ecm240', 'wzp6_ee_eeH_Hgg_ecm240', 'wzp6_ee_mumuH_Hmumu_ecm240', 'wzp6_ee_nunuH_Hmumu_ecm240', 'wzp6_ee_eeH_Hmumu_ecm240', 'wzp6_ee_mumuH_HZa_ecm240', 'wzp6_ee_nunuH_HZa_ecm240', 'wzp6_ee_eeH_HZa_ecm240', 'wzp6_ee_mumuH_Hss_ecm240', 'wzp6_ee_nunuH_Hss_ecm240', 'wzp6_ee_eeH_Hss_ecm240', 'wzp6_ee_mumuH_Haa_ecm240', 'wzp6_ee_eeH_Haa_ecm240']}
                 }
                            
                                   
                                   
                                   
                        

              


legend = {}
legend['Signal,HZZ'] = 'Signal'
legend['llH,HZZ'] = 'l^{+}l^{-}H, HZZ'
legend['WW'] = 'WW'
legend['ZZ'] = 'ZZ'
legend['Htautau'] = 'H#rightarrow#tau#tau'
legend['HWW'] = 'H#rightarrowWW'
legend['Hbb'] = 'H#rightarrowbb'
legend['H->other decays'] = 'H#rightarrowother'
#legend['Mee'] = 'Mee'
#legend['VV'] = 'VV boson'
