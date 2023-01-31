"""
25 January 2023
Abraham Tishelman-Charny

The purpose of this python module is to run the plotting step of the FCC analysis for Z(cc)H. Started with examples in repo. 

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
print("JobName:",JobName)

import ROOT

# global parameters
intLumi        = 5e+06 #in pb-1
ana_tex        = 'e^{+}e^{-} #rightarrow Z(cc)H'
delphesVersion = '3.4.2'
energy         = 240.0
collider       = 'FCC-ee'

if(EOSoutput):
    inputDir = f"/eos/user/a/atishelm/ntuples/FCC/{JobName}/final/"
else:
    inputDir       = f'{JobName}/final/'

formats        = ['png','pdf']
yaxis          = ['lin','log']
stacksig       = ['nostack']
#stacksig       = ['stack','nostack']
outdir         = f'/eos/user/a/atishelm/www/FCC/{JobName}/plots/'

variables = [
    "N_selected_jets",
    # "jets_pt",
    # "jets_y",
    # "jets_p",
    # "jets_e",

    "N_dijet_pair",
    # "dijet_pair_all_pt",    
    "dijet_pair_all_m",
    # "dijet_pair_all_recoil_m",

    # "dijet_pair_nearZpeak_pt",    
    "dijet_pair_nearZpeak_m",
    # "dijet_pair_nearZpeak_recoil_m",

]

###Dictonnary with the analysis name as a key, and the list of selections to be plotted for this analysis. The name of the selections should be the same than in the final selection
selections = {}
selections['ZccH'] = ["sel0"]
selections['ZccH_combined'] = ["sel0"]
#selections['ZccH']   = ["sel0", "sel1", "sel2", "sel3", "sel4"]
#selections['ZccH_combined']   = ["sel0", "sel1", "sel2", "sel3", "sel4"]

extralabel = {}
extralabel['sel0'] = "Jet 0 p_{T} > 1 GeV"
#extralabel['sel1'] = "Jet 0 p_{T} > 20 GeV"
#extralabel['sel2'] = "Jet 0 p_{T} > 30 GeV"
#extralabel['sel3'] = "Jet 0 p_{T} > 40 GeV"
#extralabel['sel4'] = "Jet 0 p_{T} > 50 GeV"
#extralabel['sel1'] = "Selection: N_{Z} = 1; 80 GeV < m_{Z} < 100 GeV"

colors = {}

# exclusive 

# ZccH
# colors['ZccHbb'] = ROOT.kRed
# colors['ZccHmumu'] = ROOT.kRed+2
# colors['ZccHWW'] = ROOT.kGreen
# colors['ZccHgg'] = ROOT.kGreen+4
# colors['ZccHZa'] = ROOT.kBlue
# colors['ZccHss'] = ROOT.kBlue+2
# colors['ZccHcc'] = ROOT.kMagenta-9
# colors['ZccHmumu'] = ROOT.kMagenta+2
# colors['ZccHZZ'] = ROOT.kBlack
# colors['ZccHaa'] = ROOT.kGray
# colors['ZccHtautau'] = ROOT.kViolet

#ROOT.gStyle.SetPalette(55)

colors['ZccHbb'] = 1
colors['ZccHmumu'] = 2
colors['ZccHWW'] = 3
colors['ZccHgg'] = 4
colors['ZccHZa'] = 5
colors['ZccHss'] = 6
colors['ZccHcc'] = 7
colors['ZccHmumu'] = 8
colors['ZccHZZ'] = 9
colors['ZccHaa'] = 10
colors['ZccHtautau'] = 11

#ROOT.gStyle.SetPalette(55)

print("colors:",colors)

# # https://loading.io/color/feature/Spectral-10/
# RGBs = [
#     [158, 1, 66],
#     [213, 62, 79],
#     [244, 109, 67],
#     [253, 174, 97],
#     [254, 254, 139],
#     [230, 245, 152],
#     [171, 221, 164],
#     [102, 194, 165],
#     [50, 136, 189],
#     [94, 79, 162],
# ]

# for ci, key in enumerate(colors):
#     tmp_color = ROOT.TColor(999, 0.1, 0.2, 0.3)
#     #TColor *color = new TColor(ci, 0.1, 0.2, 0.3);
#     c1, c2, c3 = RGBs[ci]
#     tmp_color = tmp_color.SetRGB(c1, c2, c3)
#     print("tmp_color:",tmp_color)
#     colors[key] = tmp_color

print("colors:",colors)

# VV 
colors['WW'] = ROOT.kBlue+1
colors['ZZ'] = ROOT.kGreen+2

# inclusive 
colors['ZccH'] = ROOT.kGreen # signal 
colors['VV'] = ROOT.kRed # background

plots = {}
plots['ZccH'] = {

    'signal' : {
        'ZccHbb' : ['wzp6_ee_ccH_Hbb_ecm240'],
        'ZccHmumu': ['wzp6_ee_ccH_Hmumu_ecm240'],
        'ZccHWW' : ['wzp6_ee_ccH_HWW_ecm240'],
        'ZccHgg' :       ['wzp6_ee_ccH_Hgg_ecm240'],
        'ZccHZa' :       ['wzp6_ee_ccH_HZa_ecm240'],
        'ZccHss':        ['wzp6_ee_ccH_Hss_ecm240'],
        'ZccHcc' :       ['wzp6_ee_ccH_Hcc_ecm240'],
        'ZccHmumu' :       ['wzp6_ee_ccH_Hmumu_ecm240'],
        'ZccHZZ':        ['wzp6_ee_ccH_HZZ_ecm240'],	
        'ZccHtautau':        ['wzp6_ee_ccH_Htautau_ecm240'],
        'ZccHaa':        ['wzp6_ee_ccH_Haa_ecm240'],
        'ZccHbb':        ['wzp6_ee_ccH_Hbb_ecm240'],
    },

    'backgrounds' : {
        'WW':['p8_ee_WW_ecm240'],
        'ZZ':['p8_ee_ZZ_ecm240']
        }
}

plots['ZccH_combined'] = {

    'signal' : {

        'ZccH' : ['wzp6_ee_ccH_Hbb_ecm240', 
                  'wzp6_ee_ccH_Hmumu_ecm240',
                  'wzp6_ee_ccH_HWW_ecm240',
                  'wzp6_ee_ccH_Hgg_ecm240',
                  'wzp6_ee_ccH_HZa_ecm240',
                  'wzp6_ee_ccH_Hss_ecm240',
                  'wzp6_ee_ccH_Hcc_ecm240',
                  'wzp6_ee_ccH_HZZ_ecm240',
                  'wzp6_ee_ccH_Htautau_ecm240',
                  'wzp6_ee_ccH_Haa_ecm240'],
    },

    'backgrounds' : { 'VV' : ['p8_ee_WW_ecm240', 'p8_ee_ZZ_ecm240']}

}

legend = {}
legend['ZccHbb'] = 'Z(cc)H(bb)'
legend['ZccHmumu'] = 'Z(cc)H(\mu\mu)'
legend['ZccHWW'] = 'Z(cc)H(WW)'
legend['ZccHgg'] = 'Z(cc)H(gg)'
legend['ZccHZa'] = 'Z(cc)H(Z\gamma)'
legend['ZccHss'] = 'Z(cc)H(ss)'
legend['ZccHcc'] = 'Z(cc)H(cc)'
legend['ZccHZZ'] = 'Z(cc)H(ZZ)'
legend['ZccHaa'] = 'Z(cc)H(\gamma\gamma)'
legend['ZccHtautau'] = 'Z(cc)H(\\tau\\tau)'

# VV 
legend['WW'] = 'WW'
legend['ZZ'] = 'ZZ'

# inclusive 
legend['ZccH'] = 'ZccH'
legend['VV'] = 'VV'