// Author: Claude Pruneau   12/12/16

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class CanvasConfiguration
 \ingroup WAC

 Utility class used to define the parameters of a root canvas
 */

#include "Event.hpp"

// Factory<Particle> * particleFactory = 0;
//
// extern Factory<Particle> * getParticleFactory(int nParticleMax=5000)
//{
//   if (!particleFactory)
//     {
//     particleFactory = new Factory<Particle>();
//     particleFactory->initialize(nParticleMax);
//     }
//   return particleFactory;
// }

ClassImp(Event)

  Event::Event()
  : index(0),
    eventNumber(0),
    nProjectile(0),
    nTarget(0),
    nParticleTotal(0),
    classestimator(kV0M),
    nParticles(0),
    multiplicityclass(0),
    multiplicity(0),
    impactParameter(0),
    other(0),
    inelgth0(false),
    V0AM(0),
    V0CM(0),
    CL1M(0),
    CL1EtaGapM(0),
    dNchdEta(0),
    fhNPartTot(nullptr),
    fhMultiplicity(nullptr),
    fhV0Multiplicity(nullptr),
    fhCL1Multiplicity(nullptr),
    fhCL1EtaGapMultiplicity(nullptr),
    fhV0MMultPercentile(nullptr),
    fhCL1MultPercentile(nullptr),
    fhCL1EtaGapMultPercentile(nullptr),
    ownParticleFactory(nullptr)
{
  fhNPartTot = new TH1F("EventNpart", "Event analyzed particles;number of particles;counts", 1000, -0.5, 1000 - 0.5);
  fhMultiplicity = new TH1F("EventMultiplicity", "Event multiplicity;multiplicity (%);counts", 101, -0.5, 101 - 0.5);
  fhV0Multiplicity = new TH2F("V0Multiplicity", "V0M;V0M;d#it{N}/d#eta;counts", 300, -9.5, 300 - 9.5, 250, -9.5, 250 - 9.5);
  fhCL1Multiplicity = new TH2F("CL1Multiplicity", "CL1M;CL1M;d#it{N}/d#eta;counts", 300, -9.5, 300 - 9.5, 250, -9.5, 250 - 9.5);
  fhCL1EtaGapMultiplicity = new TH2F("CL1EtaGapMultiplicity", "CL1M (excl |#eta|<0.8);CL1M;d#it{N}/d#eta;counts", 300, -9.5, 300 - 9.5, 250, -9.5, 250 - 9.5);
}

//////////////////////////
// DTOR
//////////////////////////
Event::~Event()
{
  delete fhNPartTot;
  delete fhMultiplicity;
  delete fhV0Multiplicity;
  delete fhCL1Multiplicity;
  delete fhCL1EtaGapMultiplicity;
  delete fhV0MMultPercentile;
  delete fhCL1MultPercentile;
  delete fhCL1EtaGapMultPercentile;
  /* null for the singleton Event::event; non-null only for instances created via createReconstructed() */
  delete ownParticleFactory;
}

////////////////////////////////////////////////////
// Call before restarting simulation
////////////////////////////////////////////////////
void Event::clear()
{
  index = 0;
  nParticles = 0;
  eventNumber = 0;
  nProjectile = 0;
  nTarget = 0;
  nParticleTotal = 0;
  nParticles = 0;
  multiplicityclass = 0;
  multiplicity = 0;
  impactParameter = 0;
  inelgth0 = false;
  other = 0;
  V0AM = 0;
  V0CM = 0;
  CL1M = 0;
  CL1EtaGapM = 0;
  dNchdEta = 0;
  particleFactory()->reset();
}

////////////////////////////////////////////////////
// Call before generating new event
////////////////////////////////////////////////////
void Event::reset()
{
  index++;
  nParticles = 0;
  eventNumber = 0;
  nProjectile = 0;
  nTarget = 0;
  nParticleTotal = 0;
  nParticles = 0;
  multiplicityclass = 0;
  multiplicity = 0;
  impactParameter = 0;
  other = 0;
  inelgth0 = false;
  V0AM = 0;
  V0CM = 0;
  CL1M = 0;
  CL1EtaGapM = 0;
  dNchdEta = 0;
  particleFactory()->reset();
}

void Event::setMultiplicityPercentiles(TFile* f)
{
  Info("Event::setMultiplicityPercentiles()", "From file %s", f->GetName());
  fhV0MMultPercentile = (TH1*)f->Get("V0MCentMult")->Clone();
  fhCL1MultPercentile = (TH1*)f->Get("CL1MCentMult")->Clone();
  fhCL1EtaGapMultPercentile = (TH1*)f->Get("CL1EtaGapMCentMult")->Clone();

  if (fhV0MMultPercentile == nullptr || fhCL1MultPercentile == nullptr || fhCL1EtaGapMultPercentile == nullptr) {
    Fatal("Event::setMultiplicityPercentiles()", "Percentiles histograms not correctly loaded. ABORTING!!!");
    return;
  }
}

void Event::saveHistograms(TDirectory* dir)
{
  dir->cd();
  fhNPartTot->Write();
  fhMultiplicity->Write();
  fhV0Multiplicity->Write();
  fhCL1Multiplicity->Write();
  fhCL1EtaGapMultiplicity->Write();
}

/// \brief all multiplicity information has been provided
/// settle the multiplicity and store the multiplicity class
/// if feasible
/// \param npart the number of selected particles
void Event::settleMultiplicity(int npart)
{
  multiplicityclass = 105;
  multiplicity = 0;
  if (inelgth0) {
    nParticles = npart;
    if (fhNPartTot != nullptr) {
      fhNPartTot->Fill(npart);
    }
    if (fhV0Multiplicity != nullptr) {
      fhV0Multiplicity->Fill(V0AM + V0CM, dNchdEta);
    }
    if (fhCL1Multiplicity != nullptr) {
      fhCL1Multiplicity->Fill(CL1M, dNchdEta);
    }
    if (fhCL1EtaGapMultiplicity != nullptr) {
      fhCL1EtaGapMultiplicity->Fill(CL1EtaGapM, dNchdEta);
    }
    switch (classestimator) {
      case kV0M:
        if (fhV0MMultPercentile != nullptr) {
          multiplicityclass = fhV0MMultPercentile->GetBinContent(fhV0MMultPercentile->FindFixBin(V0AM + V0CM));
          multiplicity = V0AM + V0CM;
        }
        break;
      case kCL1:
        if (fhCL1MultPercentile != nullptr) {
          multiplicityclass = fhCL1MultPercentile->GetBinContent(fhCL1MultPercentile->FindFixBin(CL1M));
          multiplicity = CL1M;
        }
        break;
      case kCL1GAP:
        if (fhCL1EtaGapMultPercentile != nullptr) {
          multiplicityclass = fhCL1EtaGapMultPercentile->GetBinContent(fhCL1EtaGapMultPercentile->FindFixBin(CL1EtaGapM));
          multiplicity = CL1EtaGapM;
        }
        break;
      default:
        break;
    }
    fhMultiplicity->Fill(multiplicityclass);
  }
}

///////////////////////////////////////////////////////
// Print properties of this event at the given output
///////////////////////////////////////////////////////
void Event::printProperties(ostream& output)
{
  output << "                 Event# " << index << endl;
  output << " nParticles(Accepted) : " << nParticles << endl;
  output << "          eventNumber : " << eventNumber << endl;
  output << "          nProjectile : " << nProjectile << endl;
  output << "              nTarget : " << nTarget << endl;
  output << "       nParticleTotal : " << nParticleTotal << endl;
  output << "   multiplicity class : " << multiplicityclass << endl;
  output << "         multiplicity : " << multiplicity << endl;
  output << "      impactParameter : " << impactParameter << endl;
  output << "                other : " << other << endl;
  output << "             inelgth0 : " << (inelgth0 ? "true" : "false") << endl;
  output << "                 V0AM : " << V0AM << endl;
  output << "                 V0CM : " << V0CM << endl;
  output << "                 CL1M : " << CL1M << endl;
  output << "           CL1EtaGapM : " << CL1EtaGapM << endl;
  output << "             dNchdEta : " << dNchdEta << endl;

  for (int iParticle = 0; iParticle < nParticles; iParticle++) {
    particleFactory()->getObjectAt(iParticle)->printProperties(output);
  }
}

Event* Event::event = nullptr;
Event* Event::getEvent()
{
  if (!event) {
    event = new Event();
  }
  return event;
}

/// \brief Create a non-singleton Event with its own Factory<Particle>.
///
/// Used by DetectorEffectsTask to hold the per-event particle set after
/// detector effects have been applied, side-by-side with the singleton
/// Event::getEvent() that holds the raw generated particles.
Event* Event::createReconstructed()
{
  Event* e = new Event(); /* protected ctor is reachable from this static member */
  e->ownParticleFactory = new Factory<Particle>();
  e->ownParticleFactory->initialize(Particle::factorySize);
  return e;
}

/// \brief Reserve and return the next slot in this Event's particle factory.
///
/// Does NOT increment nParticles -- the caller (DetectorEffectsTask) decides
/// how many slots are actually used and then calls setNParticlesAccepted().
Particle* Event::appendParticle()
{
  return particleFactory()->getNextObject();
}

/// \brief Copy multiplicity / event-classification info from another Event.
///
/// Used by the detector-effects pass so the reconstructed Event reports the
/// same centrality / multiplicity class as the raw event when read by
/// EventFilter::accept() and the analyzers' fillEventWiseInfo().  The
/// per-event histograms (fhNPartTot, fhMultiplicity, ...) are NOT shared --
/// they stay zero for the reconstructed Event because settleMultiplicity()
/// is not called on it.
void Event::copyEventLevelInfoFrom(const Event& src)
{
  classestimator = src.classestimator;
  multiplicityclass = src.multiplicityclass;
  multiplicity = src.multiplicity;
  impactParameter = src.impactParameter;
  other = src.other;
  inelgth0 = src.inelgth0;
  V0AM = src.V0AM;
  V0CM = src.V0CM;
  CL1M = src.CL1M;
  CL1EtaGapM = src.CL1EtaGapM;
  dNchdEta = src.dNchdEta;
  eventNumber = src.eventNumber;
  nProjectile = src.nProjectile;
  nTarget = src.nTarget;
  nParticleTotal = src.nParticleTotal;
}
