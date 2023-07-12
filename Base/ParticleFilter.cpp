//
//  ParticleFilter.cpp
//  MyMC
//
//  Created by Claude Pruneau on 12/12/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include <TMath.h>
#include "ParticleFilter.hpp"

//////////////////////////////////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
ParticleFilter<r>::ParticleFilter(ParticleFilter::SpeciesSelection pid,
                                  ParticleFilter::ChargeSelection charge,
                                  ParticleFilter::FeedDownRejection feeddown,
                                  double minPt, double maxPt,
                                  double minRapPseudo, double maxRapPseudo)
  : pidRequested(pid), chargeRequested(charge), feedDownRejection(feeddown), min_pt(minPt), max_pt(maxPt), min_rappseudo(minRapPseudo), max_rappseudo(maxRapPseudo)
{
  // no ops
}

//////////////////////////////////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
ParticleFilter<r>::~ParticleFilter()
{
  // no ops
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
TString ParticleFilter<r>::getName()
{
  TString name;
  switch (pidRequested) {
    default:
      name = "WakyWaky";
      break;
    case AllSpecies:
      name = "All";
      break;
    case Photon:
      name = "Gam";
      break;
    case Lepton:
      name = "Lep";
      break;
    case Electron:
      name = "El";
      break;
    case Muon:
      name = "Mu";
      break;
    case Hadron:
      name = "H";
      break;
    case Pion:
      name = "Pi";
      break;
    case Kaon:
      name = "Ka";
      break;
    case KaonL:
      name = "KaL";
      break;
    case Baryon:
      name = "B";
      break;
    case Proton:
      name = "Pr";
      break;
    case Lambda:
      name = "La";
      break;
    case ALambda:
      name = "ALa";
      break;
  }

  if (pidRequested == Photon)
    return name;

  switch (chargeRequested) {
    default:
      name = "WakyWakyWakyWaky";
      break;
    case AllCharges:
      name += "All";
      break; // all charges
    case Negative:
      name += "M";
      break; // minus
    case Positive:
      name += "P";
      break;
    case Charged:
      name += "C";
      break;
    case Neutral:
      name += "0";
      break;
  }
  /* TODO: for the time being we do not incorporate the feed down suppression in the name */
  return name;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
TString ParticleFilter<r>::getTitle()
{
  TString name;
  switch (pidRequested) {
    case AllSpecies:
      name = "All Specie";
      break;
    case Photon:
      name = "Gamma";
      break;
    case Lepton:
      name = "Lepton";
      break;
    case Electron:
      name = "Electron";
      break;
    case Muon:
      name = "Muon";
      break;
    case Hadron:
      name = "Hadron";
      break;
    case Pion:
      name = "Pion";
      break;
    case Kaon:
      name = "Kaon";
      break;
    case KaonL:
      name = "KaonL";
      break;
    case Baryon:
      name = "Baryon";
      break;
    case Proton:
      name = "Proton";
      break;
    case Lambda:
      name = "Lambda";
      break;
    case ALambda:
      name = "ALambda";
      break;
  }

  if (pidRequested == Photon)
    return name;

  switch (chargeRequested) {
    case AllCharges:
      name += " All Charge";
      break; // all charges
    case Negative:
      name += " Minus";
      break; // minus
    case Positive:
      name += " Plus";
      break;
    case Charged:
      name += " Charged";
      break;
    case Neutral:
      name += " Neutral";
      break;
  }
  /* TODO: for the time being we do not incorporate the feed down suppression in the title */
  return name;
}

template class ParticleFilter<AnalysisConfiguration::kRapidity>;
template class ParticleFilter<AnalysisConfiguration::kPseudorapidity>;

templateClassImp(ParticleFilter)
