// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef TwoPartDiffCorrelationAnalyzerME_hpp
#define TwoPartDiffCorrelationAnalyzerME_hpp

#include <TParameter.h>
#include "TFile.h"
#include "TList.h"
#include "Particle.hpp"
#include "Task.hpp"
#include "Event.hpp"
#include "EventPool.hpp"
#include "EventFilter.hpp"
#include "EventHistos.hpp"
#include "ParticleFilter.hpp"
#include "ParticlePairFilter.hpp"
#include "ParticleHistos.hpp"
#include "ParticlePairDiffHistos.hpp"
#include "ParticlePairDerivedDiffHistos.hpp"
#include "ParticlePairBalanceFunctionDiffHistos.hpp"
#include "ParticlePairCombinedDiffHistos.hpp"

#define EVENTPOOLSIZE 10
#define PARTICLESPEREVENT 2000

template <AnalysisConfiguration::RapidityPseudoRapidity r>
class TwoPartDiffCorrelationAnalyzerME : public Task
{
 public:
  //////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////
  TwoPartDiffCorrelationAnalyzerME(const TString& name,
                                   AnalysisConfiguration* configuration,
                                   Event* event,
                                   EventFilter* eventFilter,
                                   std::vector<ParticleFilter<r>*> particleFilters);
  virtual ~TwoPartDiffCorrelationAnalyzerME();
  virtual void execute();
  virtual void createHistograms();
  virtual void loadHistograms(TFile* inputFile);
  virtual void loadBaseHistograms(TDirectory* dir);
  virtual void saveHistograms(TDirectory* dir);
  virtual void addHistogramsToExtList(TList* list, bool all = false);
  virtual void scaleHistograms(double factor);
  virtual void calculateDerivedHistograms();

  //////////////////////////////////////////////////////////////
  // Data Members
  //////////////////////////////////////////////////////////////

  EventFilter* eventFilter;
  EventPool eventPool;
  std::vector<ParticleFilter<r>*> particleFilters;

  std::vector<TString> partNames;
  std::vector<int> nAccepted;

  EventHistos* event_Histos;
  std::vector<ParticleHistos*> particle_Histos;
  std::vector<std::vector<ParticlePairDerivedDiffHistos*>> pairs_Histos;
  std::vector<std::vector<ParticlePairCombinedDiffHistos*>> pairs_CIHistos;
  std::vector<std::vector<ParticlePairCombinedDiffHistos*>> pairs_CDHistos;
  std::vector<std::vector<ParticlePairBalanceFunctionDiffHistos*>> pairs_BFHistos;

  std::vector<std::vector<ParticlePairDerivedDiffHistos*>> pairs_Histos_me;
  std::vector<std::vector<ParticlePairCombinedDiffHistos*>> pairs_CIHistos_me;
  std::vector<std::vector<ParticlePairCombinedDiffHistos*>> pairs_CDHistos_me;
  std::vector<std::vector<ParticlePairBalanceFunctionDiffHistos*>> pairs_BFHistos_me;

  ClassDef(TwoPartDiffCorrelationAnalyzerME, 2)
};

#endif /* TwoPartDiffCorrelationAnalyzerME_hpp */
