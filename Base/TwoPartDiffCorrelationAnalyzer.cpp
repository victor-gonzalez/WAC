// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Two Particle Correlation Analyzer Task
 */

#include <TString.h>
#include "TwoPartDiffCorrelationAnalyzer.hpp"
#include "AnalysisConfiguration.hpp"

//////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
TwoPartDiffCorrelationAnalyzer<r>::TwoPartDiffCorrelationAnalyzer(const TString& name, AnalysisConfiguration* configuration, Event* event,
                                                                  EventFilter* ef, std::vector<ParticleFilter<r>*> particleFilters)
  : Task(name, configuration, event),
    eventFilter(ef),
    particleFilters(particleFilters),
    nAccepted(particleFilters.size(), 0),
    nAcceptedPairs{particleFilters.size(), std::vector<int>(particleFilters.size(), 0)},
    event_Histos(nullptr),
    particle_Histos{particleFilters.size(), nullptr},
    pairs_Histos{particleFilters.size(), {particleFilters.size(), nullptr}}
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::CTOR(...) Started." << endl;

  if (!eventFilter) {
    if (reportError())
      cout << "TwoPartDiffCorrelationAnalyzer::CTOR(...) eventFilter is null pointer." << endl;
    postTaskError();
    return;
  }

  for (auto pf : particleFilters) {
    if (pf == nullptr) {
      if (reportError())
        cout << "TwoPartDiffCorrelationAnalyzer::CTOR(...) particleFilters has a null pointer." << endl;
      postTaskError();
      return;
    }
  }

  TString newName = getName();
  newName += "_";
  newName += eventFilter->getName();
  setName(newName);

  for (auto pf : particleFilters) {
    partNames.push_back(pf->getName());
  }

  for (uint i = 0; i < particleFilters.size(); ++i) {
    pairs_CIHistos.push_back(std::vector<ParticlePairCombinedDiffHistos*>(particleFilters.size() - (i + 1), nullptr));
    pairs_CDHistos.push_back(std::vector<ParticlePairCombinedDiffHistos*>(particleFilters.size() - (i + 1), nullptr));
  }
  for (uint i = 0; i < uint(particleFilters.size() / 2); ++i) {
    pairs_BFHistos.push_back(std::vector<ParticlePairBalanceFunctionDiffHistos*>(uint(particleFilters.size() / 2), nullptr));
  }
}

//////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
TwoPartDiffCorrelationAnalyzer<r>::~TwoPartDiffCorrelationAnalyzer()
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::DTOR(...) Started" << endl;
  if (event_Histos != NULL)
    delete event_Histos;
  for (auto pf : particle_Histos) {
    delete pf;
  }
  for (auto& v : pairs_Histos) {
    for (auto ph : v) {
      delete ph;
    }
  }
  for (auto& v : pairs_CIHistos) {
    for (auto ph : v) {
      delete ph;
    }
  }
  for (auto& v : pairs_CDHistos) {
    for (auto ph : v) {
      delete ph;
    }
  }
  for (auto& v : pairs_BFHistos) {
    for (auto ph : v) {
      delete ph;
    }
  }
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::DTOR(...) Completed" << endl;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void TwoPartDiffCorrelationAnalyzer<r>::createHistograms()
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::initialize(...) started" << endl;
  AnalysisConfiguration* ac = (AnalysisConfiguration*)getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();

  event_Histos = new EventHistos("Event Histos", ac, debugLevel);

  for (uint i = 0; i < partNames.size(); ++i) {
    particle_Histos[i] = new ParticleHistos(partNames[i], ac, debugLevel);
  }

  if (ac->fillPairs) {
    for (uint i = 0; i < partNames.size(); ++i) {
      for (uint j = 0; j < partNames.size(); ++j) {
        pairs_Histos[i][j] = new ParticlePairDerivedDiffHistos(partNames[i] + partNames[j], ac, debugLevel);
      }
    }
    if (ac->calculateDerivedHistograms) {
      for (uint i = 0; i < partNames.size(); ++i) {
        for (uint j = 0; j < partNames.size(); ++j) {
          pairs_Histos[i][j]->createDerivedHistograms();
        }
      }
      for (uint i = 0; i < partNames.size(); ++i) {
        for (uint j = 0; j < partNames.size() - (i + 1); ++j) {
          pairs_CIHistos[i][j] = new ParticlePairCombinedDiffHistos(partNames[i] + partNames[j + i + 1] + "CI", ac, debugLevel);
          pairs_CDHistos[i][j] = new ParticlePairCombinedDiffHistos(partNames[i] + partNames[j + i + 1] + "CD", ac, debugLevel);
        }
      }
      for (uint i = 0; i < uint(particleFilters.size() / 2); ++i) {
        for (uint j = 0; j < uint(particleFilters.size() / 2); ++j) {
          pairs_BFHistos[i][j] = new ParticlePairBalanceFunctionDiffHistos(TString::Format("%.2s%.2s", partNames[2 * i].Data(), partNames[2 * j].Data()), ac, debugLevel);
        }
      }
    }
  }
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::createHistograms(...) completed" << endl;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void TwoPartDiffCorrelationAnalyzer<r>::loadHistograms(TFile* inputFile)
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::loadHistograms(...) Starting." << endl;
  /* first load the number of events as from the  cumulated parameter */
  TParameter<Long64_t>* par = (TParameter<Long64_t>*)inputFile->Get("NoOfEvents");
  eventsProcessed = par->GetVal();
  delete par;
  AnalysisConfiguration* analysisConfiguration = (AnalysisConfiguration*)getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();

  for (uint i = 0; i < partNames.size(); ++i) {
    particle_Histos[i] = new ParticleHistos(inputFile, partNames[i], analysisConfiguration, debugLevel);
  }
  if (analysisConfiguration->fillPairs) {
    for (uint i = 0; i < partNames.size(); ++i) {
      for (uint j = 0; j < partNames.size(); ++j) {
        pairs_Histos[i][j] = new ParticlePairDerivedDiffHistos(inputFile, partNames[i] + partNames[j], analysisConfiguration, debugLevel);
      }
    }
    if (analysisConfiguration->calculateDerivedHistograms) {
      for (uint i = 0; i < partNames.size(); ++i) {
        for (uint j = 0; j < partNames.size() - (i + 1); ++j) {
          pairs_CIHistos[i][j] = new ParticlePairCombinedDiffHistos(inputFile, partNames[i] + partNames[j + i + 1] + "CI", analysisConfiguration, debugLevel);
          pairs_CDHistos[i][j] = new ParticlePairCombinedDiffHistos(inputFile, partNames[i] + partNames[j + i + 1] + "CD", analysisConfiguration, debugLevel);
        }
      }
      for (uint i = 0; i < uint(particleFilters.size() / 2); ++i) {
        for (uint j = 0; j < uint(particleFilters.size() / 2); ++j) {
          pairs_BFHistos[i][j] = new ParticlePairBalanceFunctionDiffHistos(inputFile, TString::Format("%.2s%.2s", partNames[2 * i].Data(), partNames[2 * j].Data()), analysisConfiguration, debugLevel);
        }
      }
    }
  }
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::loadHistograms(...) Completed." << endl;
}

//////////////////////////////////////////////////////////////
// load the base histograms from given file
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void TwoPartDiffCorrelationAnalyzer<r>::loadBaseHistograms(TDirectory* dir)
{
  dir->cd();
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::loadHistograms(...) Starting." << endl;

  /* first load the number of events as from the  cumulated parameter */
  TParameter<Long64_t>* par = (TParameter<Long64_t>*)dir->Get("NoOfEvents");
  eventsProcessed = par->GetVal();
  delete par;
  AnalysisConfiguration* analysisConfiguration = (AnalysisConfiguration*)getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();

  for (uint i = 0; i < partNames.size(); ++i) {
    particle_Histos[i] = new ParticleHistos(dir, partNames[i], analysisConfiguration, debugLevel);
  }
  if (analysisConfiguration->fillPairs) {
    for (uint i = 0; i < partNames.size(); ++i) {
      for (uint j = 0; j < partNames.size(); ++j) {
        pairs_Histos[i][j] = new ParticlePairDerivedDiffHistos(dir, partNames[i] + partNames[j], analysisConfiguration, debugLevel);
      }
    }
    if (analysisConfiguration->calculateDerivedHistograms) {
      for (uint i = 0; i < partNames.size(); ++i) {
        for (uint j = 0; j < partNames.size(); ++j) {
          pairs_Histos[i][j]->createDerivedHistograms();
        }
      }
      for (uint i = 0; i < partNames.size(); ++i) {
        for (uint j = 0; j < partNames.size() - (i + 1); ++j) {
          pairs_CIHistos[i][j] = new ParticlePairCombinedDiffHistos(partNames[i] + partNames[j + i + 1] + "CI", analysisConfiguration, debugLevel);
          pairs_CDHistos[i][j] = new ParticlePairCombinedDiffHistos(partNames[i] + partNames[j + i + 1] + "CD", analysisConfiguration, debugLevel);
        }
      }
      for (uint i = 0; i < uint(particleFilters.size() / 2); ++i) {
        for (uint j = 0; j < uint(particleFilters.size() / 2); ++j) {
          pairs_BFHistos[i][j] = new ParticlePairBalanceFunctionDiffHistos(TString::Format("%.2s%.2s", partNames[2 * i].Data(), partNames[2 * j].Data()), analysisConfiguration, debugLevel);
        }
      }
    }
  }
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::loadHistograms(...) Completed." << endl;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void TwoPartDiffCorrelationAnalyzer<r>::saveHistograms(TDirectory* dir)
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) Saving Event histograms to file." << endl;
  if (!dir) {
    if (reportError())
      cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) output directory is a null  pointer." << endl;
    postTaskError();
    return;
  }
  dir->cd();

  /* first save the number of events as a cumulated parameter */
  TParameter<Long64_t>("NoOfEvents", eventsProcessed, '+').Write();
  AnalysisConfiguration* analysisConfiguration = (AnalysisConfiguration*)getTaskConfiguration();
  if (!analysisConfiguration) {
    if (reportError())
      cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) analysisConfiguration is a null  pointer." << endl;
    postTaskError();
    return;
  }

  /* now save the event histograms */
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) saving event histograms." << endl;
  event->saveHistograms(dir);

  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) saving singles." << endl;

  for (uint i = 0; i < partNames.size(); ++i) {
    particle_Histos[i]->saveHistograms(dir);
  }

  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) saving pairs -- perhaps." << endl;

  if (analysisConfiguration->fillPairs) {
    if (reportDebug())
      cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) saving calculated histograms." << endl;

    for (uint i = 0; i < partNames.size(); ++i) {
      for (uint j = 0; j < partNames.size(); ++j) {
        pairs_Histos[i][j]->saveHistograms(dir);
      }
    }
    if (analysisConfiguration->calculateDerivedHistograms) {
      for (uint i = 0; i < partNames.size(); ++i) {
        for (uint j = 0; j < partNames.size() - (i + 1); ++j) {
          pairs_CIHistos[i][j]->saveHistograms(dir);
          pairs_CDHistos[i][j]->saveHistograms(dir);
        }
      }
      for (uint i = 0; i < uint(particleFilters.size() / 2); ++i) {
        for (uint j = 0; j < uint(particleFilters.size() / 2); ++j) {
          pairs_BFHistos[i][j]->saveHistograms(dir);
        }
      }
    }
  }
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) Completed." << endl;
}

//////////////////////////////////////////////////////////////
// add histograms to an external list
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void TwoPartDiffCorrelationAnalyzer<r>::addHistogramsToExtList(TList* list, bool all)
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::addHistogramsToExtList(...) Saving Event histograms to file." << endl;

  /* first add the number of events as a cumulated parameter */
  list->Add(new TParameter<Long64_t>("NoOfEvents", eventsProcessed, '+'));
  AnalysisConfiguration* analysisConfiguration = (AnalysisConfiguration*)getTaskConfiguration();

  for (uint i = 0; i < partNames.size(); ++i) {
    particle_Histos[i]->addHistogramsToExtList(list, all);
  }

  if (analysisConfiguration->fillPairs) {
    for (uint i = 0; i < partNames.size(); ++i) {
      for (uint j = 0; j < partNames.size(); ++j) {
        pairs_Histos[i][j]->addHistogramsToExtList(list, all);
      }
    }
    if (analysisConfiguration->calculateDerivedHistograms) {
      for (uint i = 0; i < partNames.size(); ++i) {
        for (uint j = 0; j < partNames.size() - (i + 1); ++j) {
          pairs_CIHistos[i][j]->addHistogramsToExtList(list, all);
          pairs_CDHistos[i][j]->addHistogramsToExtList(list, all);
        }
      }
      for (uint i = 0; i < uint(particleFilters.size() / 2); ++i) {
        for (uint j = 0; j < uint(particleFilters.size() / 2); ++j) {
          pairs_BFHistos[i][j]->addHistogramsToExtList(list, all);
        }
      }
    }
  }
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::addHistogramsToExtList(...) Completed." << endl;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void TwoPartDiffCorrelationAnalyzer<r>::execute()
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::analyze(...) Starting" << endl;
  if (event != NULL) {
    if (reportDebug())
      cout << "TwoPartDiffCorrelationAnalyzer::analyze(...) analyzing " << event->getNParticles() << " particles" << endl;
  } else {
    if (reportError())
      cout << "TwoPartDiffCorrelationAnalyzer::analyze(...) event pointer is NULL. Abort." << endl;
    postTaskError();
    return;
  }

  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::analyze(...) check if event is acceptable" << endl;
  if (!eventFilter->accept(*event))
    return;
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::analyze(...) acceptable event" << endl;

  AnalysisConfiguration* analysisConfiguration = (AnalysisConfiguration*)getTaskConfiguration();

  if (!analysisConfiguration) {
    if (reportError())
      cout << "TwoPartDiffCorrelationAnalyzer::analyze(...) analysisConfiguration null pointer" << endl;
    postTaskError();
    return;
  }

  /* process the singles building the particle indexes to hurry up the pairs process */
  for (uint i = 0; i < partNames.size(); ++i) {
    nAccepted[i] = 0;
    for (uint j = 0; j < partNames.size(); ++j) {
      nAcceptedPairs[i][j] = 0;
    }
  }
  for (int iParticle = 0; iParticle < event->getNParticles(); iParticle++) {
    if (reportDebug())
      cout << "TwoPartDiffCorrelationAnalyzer::analyze(...) particle: " << iParticle << endl;
    Particle* particle = event->getParticleAt(iParticle);
    if (reportDebug())
      particle->printProperties(cout);
    int ixID = particle->ixID = ParticleFilter<r>::getAcceptedIndex(particleFilters, *particle);
    if (ixID < 0) {
      if (reportDebug()) {
        cout << "  rejected! " << endl;
      }
      continue;
    }
    particle_Histos[ixID]->fill<r>(*particle, 1.0);
    nAccepted[ixID] += 1;
    if (reportDebug()) {
      cout << "  accepted as: " << particleFilters[particle->ixID]->getName() << endl;
    }
  }

  /* now process pairs if required */
  if (analysisConfiguration->fillPairs) {
    for (int iParticle1 = 0; iParticle1 < event->getNParticles(); iParticle1++) {
      Particle& particle1 = *event->getParticleAt(iParticle1);
      int ixID1 = particle1.ixID;
      if (ixID1 < 0)
        continue;
      for (int iParticle2 = 0; iParticle2 < event->getNParticles(); iParticle2++) {
        if (iParticle1 == iParticle2)
          continue;
        Particle& particle2 = *event->getParticleAt(iParticle2);
        int ixID2 = particle2.ixID;
        if (ixID2 < 0)
          continue;

        pairs_Histos[ixID1][ixID2]->fill<r>(particle1, particle2, 1.0, 1.0);
        nAcceptedPairs[ixID1][ixID2] += 1;
      }
    }
  }
  for (uint i = 0; i < partNames.size(); ++i) {
    particle_Histos[i]->fillEventWiseInfo(event->getMultiplicityClass(), nAccepted[i], 1.0);
    if (analysisConfiguration->fillPairs) {
      for (uint j = 0; j < partNames.size(); ++j) {
        pairs_Histos[i][j]->fillEventWiseInfo(event->getMultiplicityClass(),nAcceptedPairs[i][j], 1.0);
      }
    }
  }
  eventsProcessed++;
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::execute() Completed" << endl;
}

//////////////////////////////////////////////////////////////
// calculate Derived Histograms
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void TwoPartDiffCorrelationAnalyzer<r>::calculateDerivedHistograms()
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::calculateDerivedHistograms() Starting" << endl;
  AnalysisConfiguration* analysisConfiguration = (AnalysisConfiguration*)getTaskConfiguration();
  for (uint i = 0; i < partNames.size(); ++i) {
    particle_Histos[i]->completeFill();
  }
  if (analysisConfiguration->fillPairs) {
    for (uint i = 0; i < partNames.size(); ++i) {
      particle_Histos[i]->calculateAverages();
    }
    for (uint i = 0; i < partNames.size(); ++i) {
      for (uint j = 0; j < partNames.size(); ++j) {
        /* TODO: we need to think about the bin width correction */
        pairs_Histos[i][j]->calculateDerivedHistograms(particle_Histos[i], particle_Histos[j], 1.0);
      }
    }
    for (uint i = 0; i < partNames.size(); ++i) {
      for (uint j = 0; j < partNames.size() - (i + 1); ++j) {
        pairs_CIHistos[i][j]->calculate(pairs_Histos[i][i], pairs_Histos[j + i + 1][j + i + 1], pairs_Histos[i][j + i + 1], pairs_Histos[j + i + 1][i], 0.25, 0.25, 0.25, 0.25);
        pairs_CDHistos[i][j]->calculate(pairs_Histos[i][i], pairs_Histos[j + i + 1][j + i + 1], pairs_Histos[i][j + i + 1], pairs_Histos[j + i + 1][i], -0.25, -0.25, 0.25, 0.25);
      }
    }
    /* let's infer how many balance functions              */
    /* only for charged so we assume half of the particles */
    /* getting rid of the fractional part                  */
    /* probably a better way would be to go through the    */
    /* filters and actually do only what is needed         */
    for (uint i = 0; i < uint(partNames.size() / 2); ++i) {
      for (uint j = 0; j < uint(partNames.size() / 2); ++j) {
        pairs_BFHistos[i][j]->calculate(pairs_Histos[2 * i][2 * j + 1], pairs_Histos[2 * i][2 * j], pairs_Histos[2 * i + 1][2 * j], pairs_Histos[2 * i + 1][2 * j + 1]);
      }
    }
  } else {
    if (reportDebug())
      cout << "TwoPartDiffCorrelationAnalyzer::calculateDerivedHistograms() Nothing to do, fill pairs not ordered" << endl;
  }
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::calculateDerivedHistograms() Completed" << endl;
}

//////////////////////////////////////////////////////////////
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled
// They are not by construction
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void TwoPartDiffCorrelationAnalyzer<r>::scaleHistograms(double factor)
{
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::scaleHistograms(..) Scale all primary histograms by " << factor << endl;
  for (uint i = 0; i < partNames.size(); ++i) {
    particle_Histos[i]->scale(factor);
  }
  AnalysisConfiguration* analysisConfiguration = (AnalysisConfiguration*)getTaskConfiguration();
  if (analysisConfiguration->fillPairs) {
    for (uint i = 0; i < partNames.size(); ++i) {
      for (uint j = 0; j < partNames.size(); ++j) {
        pairs_Histos[i][j]->scale(factor);
      }
    }
  }
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::scale(..) Completed" << endl;
}

template class TwoPartDiffCorrelationAnalyzer<AnalysisConfiguration::kRapidity>;
template class TwoPartDiffCorrelationAnalyzer<AnalysisConfiguration::kPseudorapidity>;

templateClassImp(TwoPartDiffCorrelationAnalyzer)
