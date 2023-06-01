//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBufferJSON.h>

#include "Event.hpp"
#include "AnalysisConfiguration.hpp"
#include "TwoPartDiffCorrelationAnalyzer.hpp"
#include "TwoPartDiffCorrelationAnalyzerME.hpp"
#include "ParticleAnalyzer.hpp"
#include "EventLoop.hpp"
#include "EventFilter.hpp"
#include "ParticleFilter.hpp"
#include "PythiaConfiguration.hpp"
#include "PythiaEventGenerator.hpp"
#include "PythiaAnalysisConfiguration.hpp"

int main(int argc, char* argv[])
{
  if (argc > 3 or argc < 2) {
    Fatal("main", "Wrong number of arguments. Use RunPythiaSimulationTwoParticleDiff jobix seed");
  }
  int jobix = stoi(argv[1]);
  long seed = stol(argv[2]);

  time_t begin, end; // time_t is a datatype to store time values.
  time(&begin);      // note time before execution

  // =========================
  // Short configuration from json file
  // =========================
  PythiaAnalysisConfiguration* conf = nullptr;
  std::ifstream inf("configuration.json");
  if (inf.is_open()) {
    TString json;
    std::string line;
    while (getline(inf, line)) {
      json += line;
    }
    inf.close();
    /* let's produce the configuration object out of the string */
    TBufferJSON::FromJSON(conf, json);
    /* show it to chek it */
    TString test_json = TBufferJSON::ToJSON(conf);
    std::cout << test_json << std::endl;
  } else {
    printf("ERROR: Configuration file not found. ABORTING!!!\n");
    return 0;
  }

  std::vector<float> abs_y = conf->abs_y;
  float min_pt = conf->min_pt;
  float max_pt = conf->max_pt;
  // int nBins_pt = int((max_pt - min_pt) / 0.1);

  long nEventsRequested = conf->nEventsRequested;
  int nEventsReport = conf->nEventsReport;
  MessageLogger::LogLevel repLevel = MessageLogger::Error;
  if (conf->logLevel == "info") {
    repLevel = MessageLogger::Info;
  } else if (conf->logLevel == "error") {
    repLevel = MessageLogger::Error;
  } else if (conf->logLevel == "warning") {
    repLevel = MessageLogger::Warning;
  } else if (conf->logLevel == "debug") {
    repLevel = MessageLogger::Debug;
  } else {
    Error("main", "Not proper log level. Please, fix it!!");
    return 0;
  }

  // ==========================
  // Event Section
  // ==========================
  Event* event = Event::getEvent();
  if (!conf->inputfile.empty()) {
    TFile* f = new TFile(conf->inputfile.c_str());
    if (f != nullptr && f->IsOpen()) {
      bool status = TH1::AddDirectoryStatus();
      TH1::AddDirectory(false);
      event->setMultiplicityPercentiles(f);
      f->Close();
      delete f;
      TH1::AddDirectory(status);
    } else {
      printf("ERROR: Input file %s configured but not found", conf->inputfile.c_str());
      if (f != nullptr) {
        delete f;
      }
      return 0;
    }
  }

  // ==========================
  // Generator Section
  // ==========================
  int nOptions = 0;
  TString** pythiaOptions = new TString*[50];
  for (auto& opt : conf->pythiaOptions) {
    pythiaOptions[nOptions++] = new TString(opt);
  }
  /* add the seed */
  pythiaOptions[nOptions++] = new TString("Random:setSeed = on");
  pythiaOptions[nOptions++] = new TString(TString::Format("Random:seed = %ld", seed));

  PythiaConfiguration* pc = new PythiaConfiguration(conf->projectileA,
                                                    conf->projectileB,
                                                    conf->energy, /* energy in GeV */
                                                    nOptions,
                                                    pythiaOptions);

  /* event selection at the generator level */
  EventFilter::EventSelection eventSelectionGen = EventFilter::MinBias;
  if (conf->geventfilter == "MB") {
    eventSelectionGen = EventFilter::MinBias;
  } else {
    Error("main", "Generator event selection %s still not supported by the analysis. Please, fix it!!", conf->geventfilter.c_str());
    return 0;
  }
  EventFilter* eventFilterGen = new EventFilter(eventSelectionGen, 0.0, 0.0);

  Task* generator;
  /* particle selection at the generator level */
  if (conf->gparticlefilter == "All" && conf->gchargefilter == "All") {
    if (conf->inrapidity) {
      ParticleFilter<AnalysisConfiguration::kRapidity>* particleFilterGen = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::AllSpecies,
                                                                                                                                 ParticleFilter<AnalysisConfiguration::kRapidity>::AllCharges,
                                                                                                                                 min_pt, max_pt,
                                                                                                                                 -abs_y[0], abs_y[0]);
      generator = new PythiaEventGenerator<AnalysisConfiguration::kRapidity>("PYTHIA", pc, event, eventFilterGen, particleFilterGen);
      generator->reportLevel = repLevel;
    } else {
      Error("main", "Generator in pseudorapidity still not supported by the analysis. Please, fix it!!");
      return 0;
    }
  } else {
    Error("main", "Launcher still not prepared for configuring different particles generation. Please, fix it!!");
    return 0;
  }

  // ==========================
  // Analysis Section
  // ==========================
  bool oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);

  int nAnalysisTasks = 200;
  Task** analysisTasks = new Task*[nAnalysisTasks];
  int iTask = 0;

  for (float y : abs_y) {
    double min_y = -y;
    double max_y = y;
    int nBins_y = int((max_y - min_y) / 0.1);

    AnalysisConfiguration* ac = new AnalysisConfiguration("PYTHIA", "PYTHIA", "1.0");
    ac->loadHistograms = false;
    ac->createHistograms = true;
    ac->scaleHistograms = false;
    ac->calculateDerivedHistograms = false;
    ac->saveHistograms = true;
    ac->resetHistograms = false;
    ac->clearHistograms = false;
    ac->forceHistogramsRewrite = false;
    ac->inputPath = "Input/";
    ac->rootInputFileName = "";
    ac->outputPath = "Output/";
    ac->rootOuputFileName = TString::Format("%s_%03d", TString::Format(conf->outputfname.c_str(), int(y * 10)).Data(), jobix).Data();
    ac->outputDirectory = TString::Format("%s", TString::Format(conf->outputfname.c_str(), int(y * 10)).Data()).Data();
    ac->histoBaseName = "TEST";

    ac->bin_edges_pt = conf->ptbins;
    ac->nBins_pt = conf->n_ptbins;
    ac->min_pt = min_pt;
    ac->max_pt = max_pt;
    ac->nBins_eta = nBins_y;
    ac->min_eta = min_y;
    ac->max_eta = max_y;
    ac->nBins_y = nBins_y;
    ac->min_y = min_y;
    ac->max_y = max_y;
    ac->nBins_phi = 72;
    ac->min_phi = 0.0;
    ac->max_phi = kTWOPI;

    ac->fillPairs = true;
    ac->fill3D = false;
    ac->fillPratt = false;
    if (conf->inrapidity) {
      ac->fillYorEta = AnalysisConfiguration::kRapidity;
    } else {
      Error("main", "Tasks in pseudorapidity still not supported by the analysis. Please, fix it!!");
      return 0;
    }

    /* the pairs taskname */
    TString taskName = TString::Format(conf->taskname.c_str(), "Pairs");

    /* event selection at the analysis task level */
    for (auto& ef : conf->teventfilter) {
      EventFilter* eventFilter = nullptr;
      if (ef == "MB") {
        eventFilter = new EventFilter(EventFilter::MinBias, 0, 0);
      } else {
        float min = 0;
        float max = 0;
        sscanf(ef.c_str(), "%f-%f", &min, &max);
        eventFilter = new EventFilter(EventFilter::Centrality, min, max);
      }
      /* for having the balance function correctly extracted the particle filters have to follow certain order */
      /* - charged particle should come always first                                                           */
      /* - particles of the same species, the positive has to come first and the negative immediately after    */
      /* - after the charged particles the might come any number of neutral                                    */
      /* - the balance function produced for neutrals will not have any sense                                  */
      std::vector<ParticleFilter<AnalysisConfiguration::kRapidity>*> particleFilters;
      if (conf->inrapidity) {
        for (auto& part : conf->tpairs) {
          auto filter = PythiaAnalysisConfiguration::particleFilter<AnalysisConfiguration::kRapidity>(part, ac);
          if (filter != nullptr) {
            particleFilters.push_back(filter);
          } else {
            return 0;
          }
        }
      } else {
        Error("main", "Tasks in pseudorapidity still not supported by the analysis. Please, fix it!!");
        return 0;
      }

      /* the two-particle analyzer */
      if (conf->inrapidity) {
        analysisTasks[iTask++] = new TwoPartDiffCorrelationAnalyzerME<AnalysisConfiguration::kRapidity>(taskName, ac, event, eventFilter, particleFilters);
      } else {
        Error("main", "Tasks in pseudorapidity still not supported by the analysis. Please, fix it!!");
        return 0;
      }
      analysisTasks[iTask - 1]->reportLevel = repLevel;

      /* single particle analysis filters and task if any */
      if (conf->tsingles.size() > 0) {
        int nParticleFilters = 0;
        if (conf->inrapidity) {
          TString singlesTtaskName = TString::Format(conf->taskname.c_str(), "Singles", int(y * 10)).Data();
          ParticleFilter<AnalysisConfiguration::kRapidity>** singleParticleFilters = new ParticleFilter<AnalysisConfiguration::kRapidity>*[50];
          for (auto& part : conf->tsingles) {
            auto filter = PythiaAnalysisConfiguration::particleFilter<AnalysisConfiguration::kRapidity>(part, ac);
            if (filter != nullptr) {
              singleParticleFilters[nParticleFilters++] = filter;
            } else {
              return 0;
            }
          }
          analysisTasks[iTask++] = new ParticleAnalyzer<AnalysisConfiguration::kRapidity>(singlesTtaskName, ac, event, eventFilter, nParticleFilters, singleParticleFilters);
        } else {
          Error("main", "Tasks in pseudorapidity still not supported by the analysis. Please, fix it!!");
          return 0;
        }
        analysisTasks[iTask - 1]->reportLevel = repLevel;
      }
    }
  }

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

  TH1::AddDirectory(oldstatus);
  cout << "<INFO> PYTHIA Model Analysis - Pair Differential Correlations Histograms  - Completed" << endl;
  time(&end); // note time after execution
  double difference = difftime(end, begin);
  cout << "<INFO> in " << difference << " seconds" << endl;
}
