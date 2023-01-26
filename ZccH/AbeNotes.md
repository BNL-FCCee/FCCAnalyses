# Notes

The purpose of this file is to keep track of notes used while running FCCAnalyses repository.

```
fccanalysis run examples/FCCee/higgs/mH-recoil/mumu/analysis_stage1_batch.py 
fccanalysis run examples/FCCee/higgs/mH-recoil/mumu/analysis_stage2_batch.py 
fccanalysis final examples/FCCee/higgs/mH-recoil/mumu/analysis_final.py
fccanalysis plots examples/FCCee/higgs/mH-recoil/mumu/analysis_plots_test.py

```

Can we tag Z(cc)H by the Hyy or Zqq (from H-->ZZ) invariant mass? This gives us 

For Z(cc)H, there are peaks (with only a 10 GeV pT cut) in mjj at 90 and 125 GeV, 125 GeV slightly higher (is scaling correct? Think so...at least BR.)
Peak at Z mass comes from ZccHZZ and ZccHZy. 
Peak at H mass comes from ZccHyy, ZccHmumu.
Can we use both windows as two signal regions? Diphoton jet peak actually has less background and similar amount of signal to Z mass peak. 

Would be useful to break each of these down, e.g. GEN requirement on di-jets comes from Z.
Use generator level information to see Zcc peak, Hyy peak, ...

## Setup

First time: 

source FCC_setup.sh
export FCCDICTSDIR="/afs/cern.ch/work/f/fccsw/public/FCCDicts" # if you want to use winter2023 json

After building:

source ./setup.sh
export FCCDICTSDIR="/afs/cern.ch/work/f/fccsw/public/FCCDicts" # if you want to use winter2023 json

## Z(cc)H 

export FCCDICTSDIR="/afs/cern.ch/work/f/fccsw/public/FCCDicts" # this
should we be including doScale? intLumi? in some of these steps where those are optional.

Process one file locally

fccanalysis run ZccH/ZccH_stage1.py --output ZccH_stage1.root --files-list /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_ccH_Hbb_ecm240/events_180562176.root --nevents 100
fccanalysis run ZccH/ZccH_stage2.py --output wzp6_ee_ccH_Hbb_ecm240.root --files-list ZccH/stage1/ZccH_stage1.root --nevents 100
fccanalysis final ZccH/ZccH_final.py 
fccanalysis plots ZccH/ZccH_plots.py 

Process all files locally 

fccanalysis run ZccH/ZccH_stage1.py 
fccanalysis run ZccH/ZccH_stage2.py
fccanalysis final ZccH/ZccH_final.py 
fccanalysis plots ZccH/ZccH_plots.py 

Process all files on HTCondor:

fccanalysis run ZccH/ZccH_stage1_batch.py 

## Notes

in /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_ccH_Hbb_ecm240/events_180562176.roo,  Jet#2.collectionID is 15, indicating this is the ReconstructedParticles collection of jets. 

Jet#2: ReconstructedParticles
Jet#3: ParticleIDs

Know this from podio_metadata->Scan("m_names")

***********************************
*    Row   * Instance *   m_names *
***********************************
*        0 *        0 * MissingET *
*        0 *        1 * MCRecoAss *
*        0 *        2 * ParticleI *
*        0 *        3 * magFieldB *
*        0 *        4 * TrackerHi *
*        0 *        5 * EFlowTrac *
*        0 *        6 * Calorimet *
*        0 *        7 *  Particle *
*        0 *        8 *    Photon *
*        0 *        9 * EFlowTrac *
*        0 *       10 *  Electron *
*        0 *       11 * EFlowPhot *
*        0 *       12 * EFlowNeut *
*        0 *       13 *       Jet *
*        0 *       14 * Reconstru *
*        0 *       15 *      Muon *
***********************************