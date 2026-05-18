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
#include <string>
#include <vector>
#include "TError.h"
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
  /// Identity-only acceptance: matches the (species, charge) pair but ignores
  /// feedDown and the kinematic window.  Used by DetectorEffectsTask to
  /// classify a particle into the species index defined by the analysis
  /// track-name vocabulary (e.g. `conf->tpairs`), INDEPENDENTLY of the
  /// per-analyzer acceptance cuts (which vary across the EventLoop's
  /// analyzers but should not affect the detector's species-level
  /// efficiency lookup).
  bool acceptIdentity(Particle& particle);
  TString getName();
  TString getTitle();

  static int getAcceptedIndex(std::vector<ParticleFilter*>, Particle& particle);

  ////////////////////////////////////////////////////////////////////////////
  // Track-name vocabulary helpers (single source of truth).
  //
  // The analysis-config classes (PythiaAnalysisConfiguration,
  // BestAnalysisConfiguration, ...) hold a `tpairs` vector of short
  // track-name strings ("AllP", "PiP", "La", "Gam", ...) that designate
  // species + charge combinations.  These helpers turn those names into
  // the (SpeciesSelection, ChargeSelection) the filter constructor
  // expects and provide an identity-only classification entry point used
  // by DetectorEffectsTask.
  ////////////////////////////////////////////////////////////////////////////
  static SpeciesSelection speciesFor(const std::string& name);
  static ChargeSelection chargeFor(const std::string& name);
  static FeedDownRejection feedDownFor(const std::string& name);

  /// Identity-only acceptance from a track-name string (combines
  /// speciesFor/chargeFor with the in-place identity check).  No
  /// ParticleFilter object is constructed.
  static bool acceptIdentity(const std::string& trackName, Particle& particle);
  /// Returns the index of the first track name in `trackNames` whose
  /// (species, charge) identity matches the particle, or -1 if none.
  static int getIndex(const std::vector<std::string>& trackNames, Particle& particle);

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

//////////////////////////////////////////////////////////////////////////////////////////
// Identity-only acceptance: matches species + charge while ignoring feedDown and the
// kinematic window.  Same first two switch blocks as accept(), but returns immediately
// after them rather than proceeding to the feedDown / pt / rapidity checks.
//////////////////////////////////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline bool ParticleFilter<r>::acceptIdentity(Particle& particle)
{
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
      return true;
    case Photon:
      return (pid == 22);
    case Lepton:
      return (pid == 11) || (pid == 12) || (pid == 13) || (pid == 14) || (pid == 15) || (pid == 16);
    case Electron:
      return (pid == 11);
    case Muon:
      return (pid == 13);
    case Hadron:
      return (pid == 211) || (pid == 321) || (pid == 2212);
    case Pion:
      return (pid == 111) || (pid == 211);
    case Kaon:
      return (pid == 321) || (pid == 311) || (pid == 310);
    case KaonL:
      return (pid == 130);
    case Baryon:
      return (pid == 2212) || (pid == 2112) || (pid == 3122);
    case Proton:
      return (pid == 2212);
    case Lambda:
      return (particle.pid == 3122);
    case ALambda:
      return (particle.pid == -3122);
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Track-name vocabulary helpers (single source of truth).  Each track-name string
// designates a (species, charge) pair; the suffix carries the charge:
//   P = positive   M = negative   C = charged   0 = neutral   A = all charges
// La and ALa are neutral by construction; Gam is the photon.  Unknown names abort with
// ::Fatal (matches the previous behaviour in PythiaAnalysisConfiguration::particleFilter).
//////////////////////////////////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline typename ParticleFilter<r>::SpeciesSelection ParticleFilter<r>::speciesFor(const std::string& name)
{
  if (name == "PiP" || name == "PiM" || name == "PiC" || name == "Pi0" || name == "PiA")
    return Pion;
  if (name == "KaP" || name == "KaM" || name == "KaC" || name == "Ka0" || name == "KaA")
    return Kaon;
  if (name == "PrP" || name == "PrM" || name == "PrC" || name == "PrA")
    return Proton;
  if (name == "La")
    return Lambda;
  if (name == "ALa")
    return ALambda;
  if (name == "Gam")
    return Photon;
  if (name == "AllP" || name == "AllM" || name == "AllC" || name == "All0" || name == "AllA")
    return AllSpecies;
  ::Fatal("ParticleFilter::speciesFor", "Particle species '%s' not supported. Please fix the analysis configuration.", name.c_str());
  return AllSpecies;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline typename ParticleFilter<r>::ChargeSelection ParticleFilter<r>::chargeFor(const std::string& name)
{
  if (name == "PiP" || name == "KaP" || name == "PrP" || name == "AllP")
    return Positive;
  if (name == "PiM" || name == "KaM" || name == "PrM" || name == "AllM")
    return Negative;
  if (name == "PiC" || name == "KaC" || name == "PrC" || name == "AllC")
    return Charged;
  if (name == "Pi0" || name == "Ka0" || name == "All0" || name == "La" || name == "ALa" || name == "Gam")
    return Neutral;
  if (name == "PiA" || name == "KaA" || name == "PrA" || name == "AllA")
    return AllCharges;
  ::Fatal("ParticleFilter::chargeFor", "Particle species '%s' not supported. Please fix the analysis configuration.", name.c_str());
  return AllCharges;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline typename ParticleFilter<r>::FeedDownRejection ParticleFilter<r>::feedDownFor(const std::string& name)
{
  if (name == "none")
    return None;
  if (name == "all")
    return AllResonances;
  ::Fatal("ParticleFilter::feedDownFor", "Resonance suppression '%s' not supported. Please fix the analysis configuration.", name.c_str());
  return None;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline bool ParticleFilter<r>::acceptIdentity(const std::string& trackName, Particle& particle)
{
  /* inline (species, charge) identity check without constructing a ParticleFilter */
  const ChargeSelection charge = chargeFor(trackName);
  bool chargeOk = false;
  switch (charge) {
    case AllCharges: chargeOk = true; break;
    case Negative:   chargeOk = (particle.charge <  0); break;
    case Positive:   chargeOk = (particle.charge >  0); break;
    case Charged:    chargeOk = (particle.charge != 0); break;
    case Neutral:    chargeOk = (particle.charge == 0); break;
  }
  if (!chargeOk)
    return false;
  const SpeciesSelection species = speciesFor(trackName);
  const double pid = TMath::Abs(particle.pid);
  switch (species) {
    case AllSpecies: return true;
    case Photon:     return (pid == 22);
    case Lepton:     return (pid == 11) || (pid == 12) || (pid == 13) || (pid == 14) || (pid == 15) || (pid == 16);
    case Electron:   return (pid == 11);
    case Muon:       return (pid == 13);
    case Hadron:     return (pid == 211) || (pid == 321) || (pid == 2212);
    case Pion:       return (pid == 111) || (pid == 211);
    case Kaon:       return (pid == 321) || (pid == 311) || (pid == 310);
    case KaonL:      return (pid == 130);
    case Baryon:     return (pid == 2212) || (pid == 2112) || (pid == 3122);
    case Proton:     return (pid == 2212);
    case Lambda:     return (particle.pid ==  3122);
    case ALambda:    return (particle.pid == -3122);
  }
  return false;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
inline int ParticleFilter<r>::getIndex(const std::vector<std::string>& trackNames, Particle& particle)
{
  for (uint i = 0; i < trackNames.size(); ++i) {
    if (acceptIdentity(trackNames[i], particle)) {
      return i;
    }
  }
  return -1;
}

#endif /* WAC_ParticleFilter */
