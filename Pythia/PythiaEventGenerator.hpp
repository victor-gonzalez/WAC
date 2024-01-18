// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_PythiaEventGenerator
#define WAC_PythiaEventGenerator
#include "TParticle.h"
#include "TClonesArray.h"
#include "TPythia8.h"
#include "Task.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"

template <AnalysisConfiguration::RapidityPseudoRapidity r>
class PythiaEventGenerator : public Task
{
 public:
  PythiaEventGenerator(const TString& name,
                       TaskConfiguration* configuration,
                       Event* event,
                       EventFilter* ef,
                       ParticleFilter<r>* pf);
  virtual ~PythiaEventGenerator();
  virtual void initialize();
  virtual void finalize();
  virtual void reset();
  virtual void execute();
  int FindPDGMotherLevel(TParticle,int); 
  int GetParticleIDDaughter(TParticle); 

  int nMax;                //  = 10000;
  TClonesArray* particles; // = new TClonesArray("TParticle", nMax);

  TPythia8* pythia8; // = new TPythia8();

  EventFilter* eventFilter;
  ParticleFilter<r>* particleFilter;

  ClassDef(PythiaEventGenerator, 0)
};

#endif /* WAC_PythiaEventGenerator */
