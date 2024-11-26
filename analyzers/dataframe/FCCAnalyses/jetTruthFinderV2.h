#ifndef FCCANALYZER_jetTruthFinderV2_H
#define FCCANALYZER_jetTruthFinderV2_H

#include "ROOT/RVec.hxx"

namespace FCCAnalyses {

Vec_i jetTruthFinderV2(ROOT::VecOps::RVec<TLorentzVector> jets_tlv, Vec_mc mc) {
    // jet truth=finder: match the gen-level partons (eventually with gluons) with the jet constituents
    // matching by mimimizing the sum of dr of the parton and all the jet constituents 

    Vec_tlv genQuarks; // Lorentz-vector of potential partons (gen-level)
    Vec_i genQuarks_pdgId; // corresponding PDG ID
    for(size_t i = 0; i < mc.size(); ++i) {
        int pdgid = abs(mc.at(i).PDG);
        if(pdgid > 6 and pdgid!=25) continue; // only quarks 
        //if(pdgid > 6 and pdgid != 21) continue; // only quarks and gluons
        TLorentzVector tlv;
        tlv.SetXYZM(mc.at(i).momentum.x,mc.at(i).momentum.y,mc.at(i).momentum.z,mc.at(i).mass);
        genQuarks.push_back(tlv);
        genQuarks_pdgId.push_back(mc.at(i).PDG);
    }

    Vec_i usedIdx;
    Vec_i result;
    Vec_d dr;
    for(size_t iJet = 0; iJet < jets_tlv.size(); ++iJet) {
        Vec_d dr;
        for(size_t iGen = 0; iGen < genQuarks.size(); ++iGen) {
            if(std::find(usedIdx.begin(), usedIdx.end(), iGen) != usedIdx.end()) {
                dr.push_back(1e99); // set infinite dr, skip
                continue;}
            if (genQuarks[iGen].Pt() < 0.2){
                dr.push_back(1e99); // set infinite dr, skip
                continue;}
            dr.push_back(jets_tlv[iJet].DeltaR(genQuarks[iGen]));
        }
        int minDrIdx = std::min_element(dr.begin(),dr.end()) - dr.begin();
        usedIdx.push_back(minDrIdx);
        result.push_back(genQuarks_pdgId[minDrIdx]);
    }
    return result;
}
}      
#endif
