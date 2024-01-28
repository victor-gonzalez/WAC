// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef TwoPartDiffCorrelationAnalyzer_hpp
#define TwoPartDiffCorrelationAnalyzer_hpp

#include <TParameter.h>
#include "TFile.h"
#include "TList.h"
#include "Task.hpp"
#include "Event.hpp"
#include "EventFilter.hpp"
#include "EventHistos.hpp"
#include "AnalysisConfiguration.hpp"
#include "ParticleFilter.hpp"
#include "ParticlePairFilter.hpp"
#include "ParticleHistos.hpp"
#include "ParticlePairDiffHistos.hpp"
#include "ParticlePairDerivedDiffHistos.hpp"
#include "ParticlePairBalanceFunctionDiffHistos.hpp"
#include "ParticlePairCombinedDiffHistos.hpp"

template <AnalysisConfiguration::RapidityPseudoRapidity r>
class TwoPartDiffCorrelationAnalyzer : public Task
{
 public:
  //////////////////////////////////////////////////////////////
  // CTOR
  //////////////////////////////////////////////////////////////
  TwoPartDiffCorrelationAnalyzer(const TString& name,
                                 AnalysisConfiguration* configuration,
                                 Event* event,
                                 EventFilter* eventFilter,
                                 std::vector<ParticleFilter<r>*> particleFilters);
  virtual ~TwoPartDiffCorrelationAnalyzer();
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
  std::vector<ParticleFilter<r>*> particleFilters;

  std::vector<TString> partNames;
  std::vector<int> nAccepted;
  std::vector<std::vector<int>>nAcceptedPairs;

  EventHistos* event_Histos;
  std::vector<ParticleHistos*> particle_Histos;
  std::vector<std::vector<ParticlePairDerivedDiffHistos*>> pairs_Histos;
  std::vector<std::vector<ParticlePairCombinedDiffHistos*>> pairs_CIHistos;
  std::vector<std::vector<ParticlePairCombinedDiffHistos*>> pairs_CDHistos;
  std::vector<std::vector<ParticlePairBalanceFunctionDiffHistos*>> pairs_BFHistos;

  ClassDef(TwoPartDiffCorrelationAnalyzer, 3)
};

#endif /* TwoPartDiffCorrelationAnalyzer_hpp */
