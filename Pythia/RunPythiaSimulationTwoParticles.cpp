//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include "Event.hpp"
#include "AnalysisConfiguration.hpp"
#include "TwoPartCorrelationAnalyzer.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"
#include "TwoPartCorrelationAnalyzer.hpp"
//#include <time.h>

int main(int argc, char* argv[])
{
  if (argc > 2) {
    Fatal("main", "Wrong number of arguments. Only job index expected");
  }
  int jobix = stoi(argv[1]);

  time_t begin, end; // time_t is a datatype to store time values.
  time(&begin);      // note time before execution

  long nEventsRequested = 1000000;
  int nEventsReport = 10000;

  // ==========================
  // Event Section
  // ==========================
  Event* event = Event::getEvent();

  // ==========================
  // Generator Section
  // ==========================
  int nOptions = 0;
  TString** pythiaOptions = new TString*[50];
  pythiaOptions[nOptions++] = new TString("Init:showChangedSettings = on");      // list changed settings
  pythiaOptions[nOptions++] = new TString("Init:showChangedParticleData = off"); // list changed particle data
  pythiaOptions[nOptions++] = new TString("Next:numberCount = 10000");           // print message every n events
  pythiaOptions[nOptions++] = new TString("Next:numberShowInfo = 1");            // print event information n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowProcess = 0");         // print process record n times
  pythiaOptions[nOptions++] = new TString("Next:numberShowEvent = 0");
  pythiaOptions[nOptions++] = new TString("SoftQCD:all = on"); // Allow total sigma = elastic/SD/DD/ND
                                                               // pythiaOptions[nOptions++] = new TString("HardQCD:all = on");
  PythiaConfiguration* pc = new PythiaConfiguration(2212 /* p */,
                                                    2212 /* p */,
                                                    14000.0, /* energy in GeV */
                                                    nOptions,
                                                    pythiaOptions);
  EventFilter* eventFilterGen = new EventFilter(EventFilter::MinBias, 0.0, 0.0);
  ParticleFilter<AnalysisConfiguration::kRapidity>* particleFilterGen = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Hadron,
                                                                                                                             ParticleFilter<AnalysisConfiguration::kRapidity>::Charged,
                                                                                                                             ParticleFilter<AnalysisConfiguration::kRapidity>::None,
                                                                                                                             0.2, 100.0,
                                                                                                                             -10.0, 10.0);
  PythiaEventGenerator<AnalysisConfiguration::kRapidity>* generator = new PythiaEventGenerator<AnalysisConfiguration::kRapidity>("PYTHIA", pc, event, eventFilterGen, particleFilterGen);

  // ==========================
  // Analysis Section
  // ==========================
  AnalysisConfiguration* ac = new AnalysisConfiguration("PYTHIA", "PYTHIA", "1.0");
  ac->loadHistograms = false;
  ac->createHistograms = true;
  ac->scaleHistograms = false;
  ac->calculateDerivedHistograms = false;
  ac->saveHistograms = true;
  ac->resetHistograms = false;
  ac->clearHistograms = false;
  ac->forceHistogramsRewrite = true;
  ac->inputPath = "Input/";
  ac->rootInputFileName = "";
  ac->outputPath = "Output/";
  ac->rootOuputFileName = TString::Format("PYTHIA_softOnHardOff_Pairs_%03d_", jobix).Data();
  ac->histoBaseName = "TEST";

  ac->nBins_pt = 28;
  ac->min_pt = 0.2;
  ac->max_pt = 3.0;
  ac->nBins_eta = 20;
  ac->min_eta = -1;
  ac->max_eta = 1;
  ac->nBins_y = 40;
  ac->min_y = -2;
  ac->max_y = 2;
  ac->nBins_phi = 72;
  ac->min_phi = 0.0;
  ac->max_phi = float(TMath::TwoPi());

  ac->fillPairs = true;
  ac->fill3D = false;
  ac->fillPratt = false;
  ac->fillYorEta = AnalysisConfiguration::kRapidity;

  TString taskName;
  int nAnalysisTasks = 20;
  Task** analysisTasks = new Task*[nAnalysisTasks];

  EventFilter* eventFilter = new EventFilter(EventFilter::MinBias, 0.0, 0.0);
  int nParticleFilters = 2;
  ParticleFilter<AnalysisConfiguration::kRapidity>** particleFilters = new ParticleFilter<AnalysisConfiguration::kRapidity>*[nParticleFilters];
  particleFilters[0] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);
  particleFilters[1] = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::Proton, ParticleFilter<AnalysisConfiguration::kRapidity>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt + 0.001, ac->max_pt, ac->min_y, ac->max_y);

  int iTask = 0;
  analysisTasks[iTask++] = new TwoPartCorrelationAnalyzer<AnalysisConfiguration::kRapidity>("NarrowPPPM", ac, event, eventFilter, particleFilters[0], particleFilters[1]); // P+ vs P-
  nAnalysisTasks = iTask;

  // ==========================
  // Event Loop
  // ==========================

  EventLoop* eventLoop = new EventLoop();
  eventLoop->addTask((Task*)generator);
  for (int iAnalysisTask = 0; iAnalysisTask < nAnalysisTasks; iAnalysisTask++) {
    eventLoop->addTask(analysisTasks[iAnalysisTask]);
  }
  eventLoop->run(nEventsRequested, nEventsReport);

  cout << "<INFO> PYTHIA Model Analysis - Pair Differential Correlations Histograms  - Completed" << endl;
  time(&end); // note time after execution
  double difference = difftime(end, begin);
  cout << "<INFO> in " << difference << " seconds";
}
