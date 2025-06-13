// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef BESTEVENTREADER
#define BESTEVENTREADER

#include <fstream>

#include "Task.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"

template <AnalysisConfiguration::RapidityPseudoRapidity r>
class BestEventReader : public Task
{
 public:
  BestEventReader(const TString& name,
                  TaskConfiguration* configuration,
                  Event* event,
                  EventFilter* ef,
                  ParticleFilter<r>* pf);

  virtual ~BestEventReader();
  virtual void initialize();
  virtual void finalize();
  virtual void reset();
  virtual void execute();

  ////////////////////////////////////////////////////////////////////////////////////////
  // Data members
  ////////////////////////////////////////////////////////////////////////////////////////
  ifstream fInputFile;
  int nEventsRead;

  EventFilter* eventFilter;
  ParticleFilter<r>* particleFilter;

  ClassDef(BestEventReader, 1)
};

#endif // BESTEVENTREADER
