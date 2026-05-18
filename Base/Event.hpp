// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/

#ifndef Event_hpp
#define Event_hpp
#include "Particle.hpp"
#include "TError.h"
#include "TFile.h"
#include "TH2.h"
#include "TMathBase.h"

class Event
{
 protected:
  Event();
  static Event* event;

 public:
  enum multest {
    kV0M,
    kCL1,
    kCL1GAP
  };
  virtual ~Event();
  virtual void clear();
  virtual void reset();
  virtual void printProperties(ostream& output);
  void setMultiplicityPercentiles(TFile* f);
  void saveHistograms(TDirectory* dir);

  Particle* getParticleAt(int index);
  double getMultiplicityClass() { return multiplicityclass; };
  double getMultiplicity() { return multiplicity; };
  double getNParticles() { return nParticles; }
  double getImpactParameter() { return impactParameter; }
  double getOther() { return other; }
  bool addParticleToMultiplicity(Particle& particle);
  void settleMultiplicity(int npart);

  ///////////////////////////////////////////////////////////////////////////
  // Reconstructed-event support: a non-singleton Event that owns its own
  // Factory<Particle>, used by DetectorEffectsTask to hold the per-event
  // particle set after detector effects have been applied.
  ///////////////////////////////////////////////////////////////////////////
  static Event* createReconstructed();             ///< new Event with its own particle factory
  Particle* appendParticle();                      ///< get next factory slot (does NOT bump nParticles)
  void setNParticlesAccepted(int n) { nParticles = n; }
  void copyEventLevelInfoFrom(const Event& src);   ///< copy multiplicity-related fields from src

  //////////////////////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////////////////////
 private:
  long index;
  int eventNumber;
  int nProjectile;
  int nTarget;
  int nParticleTotal;
  multest classestimator;

  int nParticles; // number of particles to analyze
  double multiplicityclass;
  double multiplicity;
  double impactParameter;
  double other;
  bool inelgth0;
  int V0AM;
  int V0CM;
  int CL1M;
  int CL1EtaGapM;
  int dNchdEta;
  TH1F* fhNPartTot;                     ///< total number of particles analyzed
  TH1F* fhMultiplicity;                 ///< the multiplicity distribution
  TH2F* fhV0Multiplicity;               ///< the V0M multiplicity histogram
  TH2F* fhCL1Multiplicity;              ///< the CL1 multiplicity histogram
  TH2F* fhCL1EtaGapMultiplicity;        ///< the CL1 with an eta gap multiplicity histogram
  const TH1* fhV0MMultPercentile;       ///< the V0M Centrality / Multiplicity percentile estimation histogram
  const TH1* fhCL1MultPercentile;       ///< the CL1 Centrality / Multiplicity percentile estimation histogram
  const TH1* fhCL1EtaGapMultPercentile; ///< the CL1 with an eta gap Centrality / Multiplicity percentile estimation histogram

  /// If non-null, this Event owns its own particle factory (used by the
  /// reconstructed Event created via createReconstructed()).  When null,
  /// the global Particle::getFactory() singleton is used (default behaviour
  /// for Event::getEvent()).
  Factory<Particle>* ownParticleFactory;

  /// Inline accessor returning the right factory for this Event instance.
  Factory<Particle>* particleFactory()
  {
    return ownParticleFactory ? ownParticleFactory : Particle::getFactory();
  }

 public:
  static Event* getEvent();

  ClassDef(Event, 1)
};

//
// Get the particle at the given index
//
inline Particle* Event::getParticleAt(int index)
{
  if (index >= 0 && index < nParticles) {
    return particleFactory()->getObjectAt(index);
  } else
    return 0;
}

/// \brief add the passed particle to the event multiplicity counting
/// \arg p the passed generated particle
/// \returns bool if the event is accepted according to the passed particle
/// false otherwise
inline bool Event::addParticleToMultiplicity(Particle& p)
{
  /* on the fly MC production */
  /* get event multiplicity according to the passed eta range  */
  /* event multiplicity as number of primary charged particles */
  /* based on AliAnalysisTaskPhiCorrelations implementation    */
  int pdgabs = TMath::Abs(p.pid);
  auto addTo = [](Particle& p, int& est, float etamin, float etamax) {
    if (p.eta < etamax && etamin < p.eta) {
      est = est + 1;
    }
  };

  /* pdg checks */
  switch (pdgabs) {
    case 9902210: /* proton diffractive */
      /* the event is rejected by setting quite negative values in the multiplicity */
      V0AM = kMinInt / 2 + 1;
      V0CM = kMinInt / 2 + 1;
      CL1M = kMinInt;
      CL1EtaGapM = kMinInt;
      dNchdEta = kMinInt;
      Error("Event::addParticleToMultiplicity", "Proton diffractive event!!! Rejecting it!!");
      return false;
      break;
    case 211:  /* pions */
    case 321:  /* kaons */
    case 2212: /* protons */
      /* not clear if we should use IsPhysicalPrimary here */
      if (0.001 < p.pt && p.pt < 50.0) {
        if (p.eta < 1.0 && -1.0 < p.eta) {
          inelgth0 = true;
        }
        addTo(p, V0AM, 2.8, 5.1);
        addTo(p, V0CM, -3.7, -1.7);
        addTo(p, CL1M, -1.4, 1.4);
        addTo(p, CL1EtaGapM, -1.4, -0.8);
        addTo(p, CL1EtaGapM, 0.8, 1.4);
        addTo(p, dNchdEta, -0.5, 0.5);
      }
      break;
    default:
      break;
  }
  return true;
}

#endif /* Event_hpp */
