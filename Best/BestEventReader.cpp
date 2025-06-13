// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#include <TError.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMath.h>
#include "BestEventReader.hpp"

template <AnalysisConfiguration::RapidityPseudoRapidity r>
BestEventReader<r>::BestEventReader(const TString& name,
                                    TaskConfiguration* configuration,
                                    Event* event,
                                    EventFilter* ef,
                                    ParticleFilter<r>* pf)
  : Task(name, configuration, event),
    eventFilter(ef),
    particleFilter(pf)
{
  if (reportDebug())
    cout << "BestEventReader::BestEventReader(...) No ops" << endl;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
BestEventReader<r>::~BestEventReader()
{
  if (reportDebug())
    cout << "BestEventReader::~BestEventReader(...) No ops" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize the event reader
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void BestEventReader<r>::initialize()
{
  if (reportDebug())
    cout << "BestEventReader::initialize() Started" << endl;

  // Read and check the file header
  enum HeaderLines {
    kFIRST,
    kSECOND,
    kTHIRD,
    kNHEADERLINES
  };

  TString currline;
  auto nextHeaderLine = [&]() {
    if (fInputFile.peek() != EOF) {
      currline.ReadLine(fInputFile);
    } else {
      ::Fatal("BestEventReader<r>::initialize()", "Premature EOF while reading header");
    }
  };
  auto checkHeaderLine = [&currline](HeaderLines nHeaderLine) {
    switch (nHeaderLine) {
      case kFIRST:
        if (currline.BeginsWith("#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge")) {
          return true;
        }
        break;
      case kSECOND:
        if (currline.BeginsWith("# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e")) {
          return true;
        }
        break;
      case kTHIRD:
        if (currline.BeginsWith("# SMASH-2.2.1")) {
          return true;
        }
        break;
      default:
        ::Fatal("BestEventReader<r>::initialize()", "Wrong header line to check. Please fix it");
        return false;
        break;
    }
    ::Fatal("BestEventReader<r>::initialize()", "Malformed header while checking line %d", nHeaderLine);
    return false;
  };

  fInputFile.open(this->taskConfiguration->dataInputFileName);
  if (!fInputFile.is_open()) {
    ::Fatal("BestEventReader<r>::initialize()", "Cannot open input file %s", taskConfiguration->dataInputFileName.Data());
  }

  for (int nLine = kFIRST; nLine < kNHEADERLINES; ++nLine) {
    nextHeaderLine();
    checkHeaderLine(static_cast<HeaderLines>(nLine));
  }

  /* the file is now left in place for reading the first and subsequent events */
  nEventsRead = 0;

  if (reportDebug())
    cout << "BestEventReader::initialize() Completed" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reset and Initialize the event reader
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void BestEventReader<r>::reset()
{
  if (reportDebug())
    cout << "BestEventReader::reset() Started" << endl;
  event->reset();
  Particle::getFactory()->reset();
  if (reportDebug())
    cout << "BestEventReader::reset() Completed" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Read a best event from file
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void BestEventReader<r>::execute()
{
  if (reportDebug())
    cout << "BestEventReader::execute() Started" << endl;

  TString currline;
  auto nextEventLine = [&]() {
    if (fInputFile.peek() != EOF) {
      currline.ReadLine(fInputFile);
      return true;
    } else {
      ::Error("BestEventReader<r>::execute()", "Premature EOF while reading event");
      return false;
    }
  };
  auto getEventInfo = [&]() {
    if (!currline.BeginsWith("# event")) {
      ::Error("execute", "Wrong start event information");
      return -1;
    } else {
      int dummyNum = 0; /* in principle should be zero */
      int numParticles = 0;
      sscanf(currline.Data(), "# event %d out %d", &dummyNum, &numParticles);
      return numParticles;
    }
  };

  Factory<Particle>* particleFactory = Particle::getFactory();
  int nParticles;
  int particlesAccepted = 0;
  int particlesCounted = 0;

  float eventPlaneAngle = TMath::TwoPi() * gRandom->Rndm();
  float cosEventPlaneAngle = std::cos(eventPlaneAngle);
  float sinEventPlaneAngle = std::sin(eventPlaneAngle);

  nextEventLine();
  nParticles = getEventInfo();

  Particle aParticle;
  for (int iParticle = 0; iParticle < nParticles; ++iParticle) {

    auto getParticle = [&currline, &sinEventPlaneAngle, &cosEventPlaneAngle](auto& aParticle) {
      enum ParticleInfo {
        kT,
        kX,
        kY,
        kZ,
        kMASS,
        kP0,
        kPX,
        kPY,
        kPZ,
        kPDG,
        kID,
        kCHARGE,
        kNOOFPARTICLEFIELDS
      };

      TObjArray* values = currline.Tokenize(" ");

      float p0 = (((TObjString*)values->At(kP0))->String()).Atof();
      float _pX = (((TObjString*)values->At(kPX))->String()).Atof();
      float _pY = (((TObjString*)values->At(kPY))->String()).Atof();
      float pX = cosEventPlaneAngle * _pX - sinEventPlaneAngle * _pY;
      float pY = sinEventPlaneAngle * _pX + cosEventPlaneAngle * _pY;
      float pZ = (((TObjString*)values->At(kPZ))->String()).Atof();
      Long64_t pdgCode = (((TObjString*)values->At(kPDG))->String()).Atoll();
      int charge = (((TObjString*)values->At(kCHARGE))->String()).Atoi();

      /* TODO: most probably we will need to incorporate the random rotatation of the event plane */
      aParticle.setPidPxPyPzE(pdgCode, charge, pX, pY, pZ, p0);

      delete values;
    };

    nextEventLine();
    getParticle(aParticle);

    /* we count the particle for multiplicity before acceptance  */
    event->addParticleToMultiplicity(aParticle);
    particlesCounted++;
    if (!particleFilter->accept(aParticle))
      continue;
    Particle* particle = particleFactory->getNextObject();
    *particle = aParticle;
    particlesAccepted++;
  }
  /* skip the final event line */
  nextEventLine();

  /* the multiplicity is settled */
  event->settleMultiplicity(particlesAccepted);
  nEventsRead++;

  if (reportDebug()) {
    cout << "BestEventReader::execute() No of accepted Particles : " << particlesAccepted << endl;
    cout << "BestEventReader::execute() No of counted Particles : " << particlesCounted << endl;
    cout << "BestEventReader::execute() Event multiplicity class : " << event->getMultiplicityClass() << endl;
    cout << "BestEventReader::execute() event completed!" << endl;
  }
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void BestEventReader<r>::finalize()
{
  if (reportDebug())
    cout << "BestEventReader::finalize() started" << endl;

  fInputFile.close();

  if (reportInfo())
    std::cout << "Read " << nEventsRead << " events" << std::endl;

  if (reportDebug())
    cout << "BestEventReader::finalize() completed" << endl;
}

template class BestEventReader<AnalysisConfiguration::kRapidity>;
template class BestEventReader<AnalysisConfiguration::kPseudorapidity>;

templateClassImp(BestEventReader);
