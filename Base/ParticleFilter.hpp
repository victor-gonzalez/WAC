// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_ParticleFilter
#define WAC_ParticleFilter
#include "TString.h"
#include "Particle.hpp"
#include "AnalysisConfiguration.hpp"

//////////////////////////////////////////////////////////////////////////////////////////
// Single Particle Filter
//
// chargeRequested:
// case -1:    accepts negative only
// case  0:    accepts neutral only
// case  1:    accepts positive only
// case  999:  accepts all
//////////////////////////////////////////////////////////////////////////////////////////

template <AnalysisConfiguration::RapidityPseudoRapidity r>
class ParticleFilter
{
 public:
  enum ChargeSelection { AllCharges,
                         Negative,
                         Positive,
                         Charged,
                         Neutral };
  enum SpeciesSelection { AllSpecies,
                          Photon,
                          Lepton,
                          Electron,
                          Muon,
                          Hadron,
                          Pion,
                          Kaon,
                          KaonL,
                          Baryon,
                          Proton,
                          Lambda,
                          ALambda };
  enum FeedDownRejection { None,
                           AllResonances };

  ParticleFilter(SpeciesSelection pidRequested,
                 ChargeSelection chargeRequested,
                 FeedDownRejection feedDownRejection,
                 double minPt,
                 double maxPt,
                 double minRapPseudo,
                 double maxRapPseudo);
  virtual ~ParticleFilter();
  bool accept(Particle& particle);
  TString getName();
  TString getTitle();

  static int getAcceptedIndex(std::vector<ParticleFilter*>, Particle& particle);

  //////////////////////////////////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////////////////////////////////

  SpeciesSelection pidRequested;
  ChargeSelection chargeRequested;
  FeedDownRejection feedDownRejection;
  double min_pt;
  double max_pt;
  double min_rappseudo;
  double max_rappseudo;

  ClassDef(ParticleFilter, 3)
};

//////////////////////////////////////////////////////////////////////////////////////////
// accept/reject the given particle based on filter parameter
//////////////////////////////////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline bool ParticleFilter<r>::accept(Particle& particle)
{
  //  enum ChargeSelection   { AllCharges, Negative, Positive, Charged, Neutral };
  //  enum SpeciesSelection  { AllSpecies, Photon, Lepton, Electron, Muon, Hadron, Pion, Kaon, Baryon, Proton, Lambda };

  bool accepting = true;
  double charge = particle.charge;
  switch (chargeRequested) {
    case AllCharges:
      accepting = true;
      break;
    case Negative:
      accepting = (charge < 0);
      break;
    case Positive:
      accepting = (charge > 0);
      break;
    case Charged:
      accepting = (charge != 0);
      break;
    case Neutral:
      accepting = (charge == 0);
      break;
  }
  if (!accepting)
    return false;
  double pid = TMath::Abs(particle.pid);
  switch (pidRequested) {
    case AllSpecies:
      accepting = true;
      break;
    case Photon:
      accepting = (pid == 22);
      break;
    case Lepton:
      accepting = (pid == 11) || (pid == 12) || (pid == 13) || (pid == 14) || (pid == 15) || (pid == 16);
      break;
    case Electron:
      accepting = (pid == 11);
      break;
    case Muon:
      accepting = (pid == 13);
      break;
    case Hadron:
      accepting = (pid == 211) || (pid == 321) || (pid == 2212);
      break;
    case Pion:
      accepting = (pid == 111) || (pid == 211);
      break;
    case Kaon:
      accepting = (pid == 321) || (pid == 311) || (pid == 310);
      break;
    case KaonL:
      accepting = (pid == 130);
      break;
    case Baryon:
      accepting = (pid == 2212) || (pid == 2112) || (pid == 3122);
      break;
    case Proton:
      accepting = (pid == 2212);
      break;
    case Lambda:
      accepting = (particle.pid == 3122);
      break;
    case ALambda:
      accepting = (particle.pid == -3122);
      break;
  }

  if (!accepting)
    return false;
  switch (feedDownRejection) {
    case None:
      accepting = true;
      break;
    case AllResonances:
      if (particle.feedDownCode != 0) {
        accepting = false;
      }
    default:
      break;
  }
  if (!accepting)
    return false;

  double pt = particle.pt;
  double eta = particle.eta;
  double y = particle.y;
  if constexpr (r == AnalysisConfiguration::kRapidity) {
    if (pt > min_pt &&
        pt <= max_pt &&
        y >= min_rappseudo &&
        y <= max_rappseudo)
      return true;
    else
      return false;
  } else {
    if (pt > min_pt &&
        pt <= max_pt &&
        eta >= min_rappseudo &&
        eta <= max_rappseudo)
      return true;
    else
      return false;
  }
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline int ParticleFilter<r>::getAcceptedIndex(std::vector<ParticleFilter*> filters, Particle& particle)
{
  for (uint i = 0; i < filters.size(); ++i) {
    if (filters[i]->accept(particle)) {
      return i;
    }
  }
  return -1;
}

#endif /* WAC_ParticleFilter */
