#ifndef  ZH_WHIZARD_ANALYZERS_H
#define  ZH_WHIZARD_ANALYZERS_H

#include <cmath>
#include <vector>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Vector3d.h"
#include "edm4hep/Vector2i.h"


/** ZH_Whizard interface.
This represents a set functions and utilities to access ZH final state, and especially classify ZHZZ events
*/
namespace FCCAnalyses{

namespace MCParticle{

  struct ZH {
    ROOT::VecOps::RVec<edm4hep::MCParticleData> H_completedecay;
    ROOT::VecOps::RVec<edm4hep::MCParticleData> Z_completedecay;
      
    ROOT::VecOps::RVec<edm4hep::MCParticleData> truth_Zq1;
    ROOT::VecOps::RVec<edm4hep::MCParticleData> truth_Zq2;
    ROOT::VecOps::RVec<edm4hep::MCParticleData> truth_Hq1;
    ROOT::VecOps::RVec<edm4hep::MCParticleData> truth_Hq2;
      
  };

  struct thetaphi {
    
    ROOT::VecOps::RVec<float> H_theta;
    ROOT::VecOps::RVec<float> H_phi;
    ROOT::VecOps::RVec<float> H_energy;
    ROOT::VecOps::RVec<float> Z_theta;
    ROOT::VecOps::RVec<float> Z_phi;
    ROOT::VecOps::RVec<float> Z_energy;
      
    ROOT::VecOps::RVec<float> truth_Zq1_theta;
    ROOT::VecOps::RVec<float> truth_Zq1_phi;
      
    ROOT::VecOps::RVec<float> truth_Zq2_theta;
    ROOT::VecOps::RVec<float> truth_Zq2_phi;
      
    ROOT::VecOps::RVec<float> truth_Hq1_theta;
    ROOT::VecOps::RVec<float> truth_Hq1_phi;
      
    ROOT::VecOps::RVec<float> truth_Hq2_theta;
    ROOT::VecOps::RVec<float> truth_Hq2_phi;
      
  };

  thetaphi fill_thetaphi_ZHdecay(const ZH &HZZ);

  // fill a struct with information from ZH, H->ZZ decays
  ZH fill_ZH_decay(const ROOT::VecOps::RVec<edm4hep::MCParticleData> &in,
                       const ROOT::VecOps::RVec<int> &ind);

  // return the invariant mass of particles in the input collection
  float invariant_mass(const ROOT::VecOps::RVec<edm4hep::MCParticleData> &in);


}//end NS MCParticle

}//end NS FCCAnalyses
#endif
