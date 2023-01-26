# Notes

The purpose of this file is to keep track of files, notes, commands used while running ZccH analysis using the FCCAnalyses repository.

## Setup 

First time: 

```
source FCC_setup.sh
```

After building:

```
source ./setup.sh
```

If one wants to use json information for a certain campaign, one should specify it with the `FCCDICTSDIR` environment variable:

```
export FCCDICTSDIR="/afs/cern.ch/work/f/fccsw/public/FCCDicts" 
```

## Z(cc)H 

Process one file locally

```
fccanalysis run ZccH/ZccH_stage1.py --output ZccH_stage1.root --files-list /eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_ccH_Hbb_ecm240/events_180562176.root --nevents 100
fccanalysis run ZccH/ZccH_stage2.py --output wzp6_ee_ccH_Hbb_ecm240.root --files-list ZccH/stage1/ZccH_stage1.root --nevents 100
fccanalysis final ZccH/ZccH_final.py 
fccanalysis plots ZccH/ZccH_plots.py 
```

To process all files locally, set `ZccH/RunConfig.yaml` variable `batch` to `0`, and run:

```
fccanalysis run ZccH/ZccH_stage1.py 
fccanalysis run ZccH/ZccH_stage2.py
fccanalysis final ZccH/ZccH_final.py 
fccanalysis plots ZccH/ZccH_plots.py 
```

Process all files on HTCondor, set `batch` to `1` and run the same commands. One can also output files to EOS, such as their CERNbox, by setting the `ZccH/RunConfig.yaml`
variable `EOSoutput` to `1`. 

## Notes

in `/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_ccH_Hbb_ecm240/events_180562176.root`,  Jet#2.collectionID is 15, indicating this is the ReconstructedParticles collection of jets. 

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

Can also tell by printing the collectionID:

```
Jet#2: ReconstructedParticles
Jet#3: ParticleIDs
```