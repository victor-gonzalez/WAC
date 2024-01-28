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

#include "ParticleAnalyzer.hpp"
#include "AnalysisConfiguration.hpp"

template <AnalysisConfiguration::RapidityPseudoRapidity r>
ParticleAnalyzer<r>::ParticleAnalyzer(const TString& name,
                                      TaskConfiguration* configuration,
                                      Event* event, EventFilter* _eventFilter,
                                      int _nParticleFilters,
                                      ParticleFilter<r>** _particleFilters)
  : Task(name, configuration, event), nParticleFilters(_nParticleFilters), eventFilter(_eventFilter), particleFilters(_particleFilters), particleHistos(nullptr), partNames(nullptr)
{
  cout << "ParticleAnalyzer::CTOR(...) Started." << endl;
  particleHistos = new ParticleHistos*[nParticleFilters];
  partNames = new TString*[nParticleFilters];
  nAccepted = new double[nParticleFilters];
  if (!eventFilter) {
    if (reportError())
      cout << "ParticleAnalyzer::CTOR(...) eventFilter is null pointer." << endl;
    postTaskError();
    return;
  }

  if (nParticleFilters < 1) {
    if (reportError())
      cout << "ParticleAnalyzer::CTOR(...) nParticleFilters<1." << endl;
    postTaskError();
    return;
  }
  if (!particleFilters) {
    if (reportError())
      cout << "ParticleAnalyzer::CTOR(...) particleFilters is null pointer." << endl;
    postTaskError();
    return;
  }

  // for (int iFilter=0; iFilter<nParticleFilters; iFilter++ )
  for (int iFilter = 0; iFilter < nParticleFilters; iFilter++) {
    if (!particleFilters[iFilter]) {
      if (reportError())
        cout << "ParticleAnalyzer::CTOR(...) particleFilter[" << iFilter << "] is a null pointer." << endl;
      postTaskError();
      return;
    } else {
      partNames[iFilter] = new TString(particleFilters[iFilter]->getName());
    }
  }
  TString newName = getName();
  newName += "_";
  newName += eventFilter->getName();
  setName(newName);
}

//////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
ParticleAnalyzer<r>::~ParticleAnalyzer()
{
  if (reportDebug())
    cout << "ParticleAnalyzer::DTOR(...) Started" << endl;
  if (particleHistos != NULL)
    delete[] particleHistos;
  if (partNames != NULL)
    delete[] partNames;
  if (reportDebug())
    cout << "ParticleAnalyzer::DTOR(...) Completed" << endl;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void ParticleAnalyzer<r>::createHistograms()
{
  if (reportDebug())
    cout << "ParticleAnalyzer::createHistograms(...) started" << endl;
  AnalysisConfiguration* ac = (AnalysisConfiguration*)getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString histoName;
  if (reportInfo())
    cout << "ParticleAnalyzer::createHistograms(...) Creating histograms for nParticleFilters:" << nParticleFilters << endl;
  for (int iFilter = 0; iFilter < nParticleFilters; iFilter++) {
    histoName = *partNames[iFilter];
    particleHistos[iFilter] = new ParticleHistos(histoName, ac, debugLevel);
  }
  if (reportDebug())
    cout << "ParticleAnalyzer::createHistograms(...) completed" << endl;
}

//////////////////////////////////////////////////////////////
// load histograms from given files
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void ParticleAnalyzer<r>::loadHistograms(TFile* inputFile)
{
  if (reportDebug())
    cout << "ParticleAnalyzer::loadHistograms(...) Starting." << endl;
  /* first load the number of events as from the  cumulated parameter */
  TParameter<Long64_t>* par = (TParameter<Long64_t>*)inputFile->Get("NoOfEvents");
  eventsProcessed = par->GetVal();
  delete par;
  AnalysisConfiguration* ac = (AnalysisConfiguration*)getTaskConfiguration();
  LogLevel debugLevel = getReportLevel();
  TString histoName;
  for (int iFilter = 0; iFilter < nParticleFilters; iFilter++) {
    histoName = *partNames[iFilter];
    particleHistos[iFilter] = new ParticleHistos(inputFile, histoName, ac, debugLevel);
  }
  if (reportDebug())
    cout << "ParticleAnalyzer::loadHistograms(...) Completed." << endl;
}

//////////////////////////////////////////////////////////////
// save histograms to given files
//////////////////////////////////////////////////////////////
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void ParticleAnalyzer<r>::saveHistograms(TDirectory* dir)
{
  if (reportDebug())
    cout << "ParticleAnalyzer::saveHistograms(...) Saving Event histograms to directory." << endl;
  if (!dir) {
    if (reportError())
      cout << "ParticleAnalyzer::saveHistograms(...) output directory is a null  pointer." << endl;
    postTaskError();
    return;
  }
  dir->cd();

  /* first save the number of events as a cumulated parameter */
  TParameter<Long64_t>("NoOfEvents", eventsProcessed, '+').Write();

  /* now save the event histograms */
  if (reportDebug())
    cout << "PartAnalyzer::saveHistograms(...) saving event histograms." << endl;
  event->saveHistograms(dir);

  for (int iFilter = 0; iFilter < nParticleFilters; iFilter++) {
    particleHistos[iFilter]->saveHistograms(dir);
  }
  if (reportDebug())
    cout << "ParticleAnalyzer::saveHistograms(...) Completed." << endl;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void ParticleAnalyzer<r>::execute()
{
  // if (reportDebug())  cout << "ParticleAnalyzer::execute(...) Starting" << endl;
  if (event != NULL) {
    if (reportDebug())
      cout << "ParticleAnalyzer::execute(...) analyzing " << event->getNParticles() << " particles" << endl;
  } else {
    if (reportError())
      cout << "ParticleAnalyzer::execute(...) event pointer is NULL. Abort." << endl;
    postTaskError();
    return;
  }

  // if (reportDebug()) cout <<"ParticleAnalyzer::analyze(...) check if event is acceptable" << endl;
  if (!eventFilter->accept(*event))
    return;
  // if (reportDebug()) cout <<"ParticleAnalyzer::analyze(...) acceptable event" << endl;
  eventsProcessed++; // count events used for fill histograms.
  bool accept;
  for (int iFilter = 0; iFilter < nParticleFilters; iFilter++)
    nAccepted[iFilter] = 0;
  for (int iParticle = 0; iParticle < event->getNParticles(); iParticle++) {
    Particle& particle = *event->getParticleAt(iParticle);
    for (int iFilter = 0; iFilter < nParticleFilters; iFilter++) {
      accept = particleFilters[iFilter]->accept(particle);
      if (accept) {
        nAccepted[iFilter]++;
        particleHistos[iFilter]->fill<r>(particle, 1.0);
      }
    }
  }
  for (int iFilter = 0; iFilter < nParticleFilters; iFilter++) {
    particleHistos[iFilter]->fillEventWiseInfo(event->getMultiplicityClass(), nAccepted[iFilter], 1.0);
  }
  // if (reportDebug()) cout << "ParticleAnalyzer::execute() Completed" << endl;
}

// =========================================================
// Scale all filled histograms by the given factor
// Derived histograms are *NOT* scaled.
// =========================================================
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void ParticleAnalyzer<r>::scaleHistograms(double factor)
{
  if (reportDebug())
    cout << "ParticleAnalyzer::scaleHistograms(..) Scale all primary histograms by " << factor << endl;
  for (int iFilter = 0; iFilter < nParticleFilters; iFilter++) {
    particleHistos[iFilter]->completeFill();
    particleHistos[iFilter]->scale(factor);
  }
  if (reportDebug())
    cout << "ParticleAnalyzer::scale(..) Completed" << endl;
}

template class ParticleAnalyzer<AnalysisConfiguration::kRapidity>;
template class ParticleAnalyzer<AnalysisConfiguration::kPseudorapidity>;

templateClassImp(ParticleAnalyzer)
