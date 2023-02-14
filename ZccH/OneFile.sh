fccanalysis run ZccH/ZccH_stage1.py --output ZccH_stage1.root --files-list /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_ccH_Hbb_ecm240/events_180562176.root --nevents 100
fccanalysis run ZccH/ZccH_stage1_Reclustering.py --output ZccH_stage1.root --files-list /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_ccH_Hbb_ecm240/events_180562176.root --nevents 100
fccanalysis run ZccH/ZccH_stage2.py --output wzp6_ee_ccH_Hbb_ecm240.root --files-list ZccH/stage1/ZccH_stage1.root --nevents 100
fccanalysis final ZccH/ZccH_final.py 
fccanalysis plots ZccH/ZccH_plots.py 
