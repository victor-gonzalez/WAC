//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <time.h>
#include <TStyle.h>
#include <TROOT.h>
#include "Event.hpp"
#include "AnalysisConfiguration.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "BestConfiguration.hpp"
#include "BestEventReader.hpp"
#include "ParticleAnalyzer.hpp"

int main(int argc, char* argv[])
{
  if (argc != 3) {
    ::Fatal("main", "Wrong number of parameters. Use RunBestSimulationSingleParticles filename nevents");
  }
  const char* filename = argv[1];
  int nEvents = stoi(argv[2]);

  time_t begin, end; // time_t is a datatype to store time values.
  time(&begin);      // note time before execution
  std::cout << "<INFO> Best reader analysis - Single Particle Histograms" << std::endl;

  long nEventsRequested = nEvents;
  int nEventsReport = std::min(nEvents, 100);

  // ==========================
  // Event Section
  // ==========================
  Event* event = Event::getEvent();

  BestConfiguration* bc = new BestConfiguration();
  bc->dataInputFileName = filename;
  EventFilter* eventFilterGen = new EventFilter(EventFilter::MinBias, 0.0, 0.0);
  ParticleFilter<AnalysisConfiguration::kRapidity>* particleFilterGen = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Hadron,
                                                                                                                             ParticleFilter<AnalysisConfiguration::kRapidity>::Charged,
                                                                                                                             ParticleFilter<AnalysisConfiguration::kRapidity>::None,
                                                                                                                             0.2, 100.0,
                                                                                                                             -10.0, 10.0);
  Task* eventReader = new BestEventReader<AnalysisConfiguration::kRapidity>("Best_afterburner", bc, event, eventFilterGen, particleFilterGen);

  // ==========================
  // Analysis Section
  // ==========================
  AnalysisConfiguration* ac = new AnalysisConfiguration("Best", "Best", "1.0");
  ac->loadHistograms = false;
  ac->createHistograms = true;
  ac->scaleHistograms = true;
  ac->calculateDerivedHistograms = false;
  ac->saveHistograms = true;
  ac->resetHistograms = false;
  ac->clearHistograms = false;
  ac->forceHistogramsRewrite = false;
  ac->inputPath = "Input/";
  ac->rootInputFileName = "";
  ac->outputPath = "Output/";
  ac->rootOuputFileName = "Best_Singles";
  ac->histoBaseName = "Base";

  ac->fillYorEta = AnalysisConfiguration::kRapidity;

  ac->nBins_pt = 100;
  ac->min_pt = 0.0;
  ac->max_pt = 100.0;
  ac->nBins_eta = 20;
  ac->min_eta = -1;
  ac->max_eta = 1;
  ac->nBins_y = 40;
  ac->min_y = -2;
  ac->max_y = 2;
  ac->nBins_phi = 36;
  ac->min_phi = 0.0;
  ac->max_phi = 2.0 * 3.1415927;

  AnalysisConfiguration* acWide = new AnalysisConfiguration(*ac);
  acWide->name = "BaseWide";
  acWide->histoBaseName = "BaseWide";
  acWide->nBins_eta = 200;
  acWide->min_eta = -10;
  acWide->max_eta = 10;
  acWide->nBins_y = 200;
  acWide->min_y = -10;
  acWide->max_y = 10;

  bool oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);

  TString taskName;
  std::vector<Task*> analysisTasks;

  EventFilter* eventFilter = new EventFilter(EventFilter::MinBias, 0.0, 0.0);
  int nParticleFilters = 12;
  ParticleFilter<AnalysisConfiguration::kRapidity>** particleFilters = new ParticleFilter<AnalysisConfiguration::kRapidity>*[nParticleFilters];
  particleFilters[0] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Hadron, ParticleFilter<AnalysisConfiguration::kRapidity>::Charged, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[1] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Hadron, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[2] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Hadron, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[3] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Pion, ParticleFilter<AnalysisConfiguration::kRapidity>::Charged, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[4] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Pion, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[5] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Pion, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[6] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Kaon, ParticleFilter<AnalysisConfiguration::kRapidity>::Charged, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[7] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Kaon, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[8] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Kaon, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[9] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Charged, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[10] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[11] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  analysisTasks.push_back(new ParticleAnalyzer<AnalysisConfiguration::kRapidity>("Narrow", ac, event, eventFilter, nParticleFilters, particleFilters));

  ParticleFilter<AnalysisConfiguration::kRapidity>** particleFiltersWide = new ParticleFilter<AnalysisConfiguration::kRapidity>*[nParticleFilters];
  particleFiltersWide[0] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Hadron, ParticleFilter<AnalysisConfiguration::kRapidity>::Charged, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[1] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Hadron, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[2] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Hadron, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[3] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Pion, ParticleFilter<AnalysisConfiguration::kRapidity>::Charged, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[4] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Pion, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[5] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Pion, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[6] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Kaon, ParticleFilter<AnalysisConfiguration::kRapidity>::Charged, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[7] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Kaon, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[8] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Kaon, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[9] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Charged, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[10] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  particleFiltersWide[11] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, acWide->min_pt + 0.001, acWide->max_pt, acWide->min_y, acWide->max_y);
  analysisTasks.push_back(new ParticleAnalyzer<AnalysisConfiguration::kRapidity>("Wide", acWide, event, eventFilter, nParticleFilters, particleFiltersWide));

  // ==========================
  // Event Loop
  // ==========================

  EventLoop* eventLoop = new EventLoop();
  eventLoop->addTask(eventReader);
  for (auto task : analysisTasks) {
    eventLoop->addTask(task);
  }
  eventLoop->run(nEventsRequested, nEventsReport);

  TH1::AddDirectory(oldstatus);

  cout << "<INFO> Best event reader analysis - Completed" << endl;
  time(&end); // note time after execution
  double difference = difftime(end, begin);
  cout << "<INFO> in " << difference << " seconds" << endl;
}
