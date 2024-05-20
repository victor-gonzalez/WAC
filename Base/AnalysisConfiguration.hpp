// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_AnalysisConfiguration
#define WAC_AnalysisConfiguration
#include "TMath.h"
#include "TaskConfiguration.hpp"
#include "Particle.hpp"
#include "EventPool.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Configuration of a given analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////
class AnalysisConfiguration : public TaskConfiguration
{
 public:
  AnalysisConfiguration(const TString& name,
                        const TString& type,
                        const TString& version);
  AnalysisConfiguration(const AnalysisConfiguration& source);
  virtual ~AnalysisConfiguration() {}

  AnalysisConfiguration& operator=(const AnalysisConfiguration& source);

  typedef enum {
    kRapidity = 0,
    kPseudorapidity = 1
  } RapidityPseudoRapidity;

  void printConfiguration(ostream& os);
  template <RapidityPseudoRapidity r>
  int getIxYEtaPhi(float eta, float phi);
  float getDphiShifted(float dphi);
  template <RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
  int getDeltaYEtaIndex(ParticleType1& particle1, ParticleType2& particle2);
  template <typename ParticleType1, typename ParticleType2>
  int getDeltaPhiIndex(ParticleType1& particle1, ParticleType2& particle2);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////
  TString histoBaseName;

  // event wise parameters
  int nBins_nPartTot;
  double min_nPartTot;
  double max_nPartTot;
  int nBins_nCell;
  double min_nCell;
  double max_nCell;
  int nBins_nPartPerCell;
  double min_nPartPerCell;
  double max_nPartPerCell;

  // single particle
  int nBins_pt;
  double min_pt;
  double max_pt;
  double range_pt;
  std::vector<double> bin_edges_pt;
  int nBins_eta;
  float min_eta;
  float max_eta;
  float range_eta;
  int nBins_y;
  float min_y;
  float max_y;
  float range_y;
  int nBins_phi;
  float min_phi;
  float max_phi;
  float range_phi;
  int nBins_phiEta;
  int nBins_phiEtaPt;
  int nBins_phiY;
  int nBins_phiYPt;

  // pair multiplicity bin correction
  float binCorrPM;
  float binCorrMP;
  float binCorrMM;
  float binCorrPP;

  // Pratt's Qinv, Qlong, Qside, Qout, additionally in kt
  int nBins_Qinv;
  double min_Qinv;
  double max_Qinv;
  int nBins_Qlong;
  double min_Qlong;
  double max_Qlong;
  int nBins_Qside;
  double min_Qside;
  double max_Qside;
  int nBins_Qout;
  double min_Qout;
  double max_Qout;
  int nBins_kT;
  double min_kT;
  double max_kT;

  int nBins_Dphi;
  float min_Dphi;
  float max_Dphi;
  float width_Dphi;
  int nBins_Dphi_shft;
  float min_Dphi_shft;
  float max_Dphi_shft;
  int nBins_Deta;
  float min_Deta;
  float max_Deta;
  float width_Deta;
  int nBins_Dy;
  float min_Dy;
  float max_Dy;
  float width_Dy;

  bool fillPairs;
  bool fillPratt;
  bool fill3D;
  RapidityPseudoRapidity fillYorEta;

  bool nuDynVsMult;
  bool nuDynVsCent;
  int nBins_mult;
  double min_mult;
  double max_mult;
  int nBins_cent;
  double min_cent;
  double max_cent;

  ClassDef(AnalysisConfiguration, 1)
};

template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline int AnalysisConfiguration::getIxYEtaPhi(float yoreta, float phi)
{

  if constexpr (r == AnalysisConfiguration::kRapidity) {
    if (!(yoreta < min_y || yoreta > max_y)) {
      int iY = int(float(nBins_y) * (yoreta - min_y) / range_y);
      int iPhi = int(float(nBins_phi) * (phi - min_phi) / range_phi);

      if (iY >= 0 && iPhi >= 0 && iY < nBins_y && iPhi < nBins_phi) {
        return nBins_phi * iY + iPhi;
      }
      return -1;
    }
    return -1;
  } else {
    if (!(yoreta < min_eta || yoreta > max_eta)) {
      int iEta = int(float(nBins_eta) * (yoreta - min_eta) / range_eta);
      int iPhi = int(float(nBins_phi) * (phi - min_phi) / range_phi);

      if (iEta >= 0 && iPhi >= 0 && iEta < nBins_eta && iPhi < nBins_phi) {
        return nBins_phi * iEta + iPhi;
      }
      return -1;
    }
  }
  return -1;
}

float const kPI = TMath::Pi();
float const kTWOPI = TMath::TwoPi();

inline float AnalysisConfiguration::getDphiShifted(float dphi)
{
  while (dphi >= max_Dphi_shft)
    dphi -= kTWOPI;
  while (dphi < min_Dphi_shft)
    dphi += kTWOPI;
  return dphi;
}

/// WARNING: for performance reasons no checks are done about the consistency
/// of particles' eta and phi within the corresponding ranges so, it is suppossed
/// the particles have been accepted and they are within that ranges
/// IF THAT IS NOT THE CASE THE ROUTINE WILL PRODUCE NONSENSE RESULTS
template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
inline int AnalysisConfiguration::getDeltaYEtaIndex(ParticleType1& particle1, ParticleType2& particle2)
{
  if constexpr (r == AnalysisConfiguration::kRapidity) {
    int ixY1 = int(float(nBins_y) * (particle1.y - min_y) / range_y);
    int ixY2 = int(float(nBins_y) * (particle2.y - min_y) / range_y);

    return ixY1 - ixY2 + nBins_y - 1;
  } else {
    int ixEta1 = int(float(nBins_eta) * (particle1.eta - min_eta) / range_eta);
    int ixEta2 = int(float(nBins_eta) * (particle2.eta - min_eta) / range_eta);

    return ixEta1 - ixEta2 + nBins_eta - 1;
  }
}

template <typename ParticleType1, typename ParticleType2>
inline int AnalysisConfiguration::getDeltaPhiIndex(ParticleType1& particle1, ParticleType2& particle2)
{
  int ixPhi1 = int(float(nBins_phi) * (particle1.phi - min_phi) / range_phi);
  int ixPhi2 = int(float(nBins_phi) * (particle2.phi - min_phi) / range_phi);

  int deltaphi_ix = ixPhi1 - ixPhi2;
  if (deltaphi_ix < 0) {
    deltaphi_ix += nBins_phi;
  }

  return deltaphi_ix;
}

#endif /* WAC_AnalysisConfiguration */
