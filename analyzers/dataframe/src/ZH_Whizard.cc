#include "FCCAnalyses/ZH_Whizard.h"
#include "FCCAnalyses/MCParticle.h"
#include <iostream>
#include <algorithm>
#include <set>


namespace FCCAnalyses{

namespace MCParticle{

  ZH fill_ZH_decay(const ROOT::VecOps::RVec<edm4hep::MCParticleData> &in,
                       const ROOT::VecOps::RVec<int> &ind) {
    ZH res;

    // Look for a Higgs boson
    ROOT::VecOps::RVec<int> Hbb_indices = get_indices(25, {5, -5}, false, false, false, false)(in, ind);
      
    if (Hbb_indices.empty()) {
      return res;
    }

    int ind_H = Hbb_indices[0];
      
    int ind_Z = 8; // assumes certain process (check MC particle visual)
    
    // Get Higgs decay products
    std::vector<int> H_idxstable = get_list_of_stable_particles_from_decay(ind_H, in, ind);
      
    std::set<int> H_idxstable_uniq(H_idxstable.begin(), H_idxstable.end());
    for (int idx : H_idxstable_uniq) {
      res.H_completedecay.push_back(in[idx]);
    }
      
    // Get Z decay products
    std::vector<int> Z_idxstable = get_list_of_stable_particles_from_decay(ind_Z, in, ind);
      
    std::set<int> Z_idxstable_uniq(Z_idxstable.begin(), Z_idxstable.end());
    for (int idx : Z_idxstable_uniq) {
      res.Z_completedecay.push_back(in[idx]);
    }
    
    return res;
  }

  
  thetaphi fill_thetaphi_ZHdecay(const ZH &HZZ){
    thetaphi res; 
    ROOT::VecOps::RVec<edm4hep::MCParticleData> H_finaldecay;
    ROOT::VecOps::RVec<edm4hep::MCParticleData> Z_finaldecay;

    // Save Higgs info
    for (auto & p: HZZ.H_completedecay) {
      if (p.generatorStatus == 1) {
	    H_finaldecay.push_back(p);
      } 
    }
 
    res.H_theta = get_theta(H_finaldecay);
    res.H_phi = get_phi(H_finaldecay);
    res.H_energy = get_e(H_finaldecay);
      
    // Save Z info
    for (auto & p: HZZ.Z_completedecay) {
      if (p.generatorStatus == 1) {
	    Z_finaldecay.push_back(p);
      } 
    }
 
    res.Z_theta = get_theta(Z_finaldecay);
    res.Z_phi = get_phi(Z_finaldecay);
    res.Z_energy = get_e(Z_finaldecay);

    return res;
  }
  

  float invariant_mass(const ROOT::VecOps::RVec<edm4hep::MCParticleData> &in) {
    TLorentzVector res;
    for (auto & p: in) {
      TLorentzVector tlv;
      tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      res += tlv;
    }
    return res.M();
  }

}//end NS MCParticle

}//end NS FCCAnalyses
