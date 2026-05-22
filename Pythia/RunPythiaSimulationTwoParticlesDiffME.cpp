//  Created by Claude Pruneau on 6/19/2020.
//  Copyright © 2020 Claude Pruneau. All rights reserved.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>
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
#include "DetectorEffectsTask.hpp"

/* generous cap: with detector effects on, each configureTasks() call appends */
/* up to 6 analyzer tasks (2 raw + 2 corrected + 2 uncorrected) instead of 2; */
/* the EventLoop's TaskCollection capacity is bumped to match in EventLoop.cpp.*/
int nAnalysisTasks = 600;
Task** analysisTasks;
int iTask = 0;

template <AnalysisConfiguration::RapidityPseudoRapidity r, AnalysisConfiguration::FillPairOptions options>
bool configureTasks(std::string efd,
                    const PythiaAnalysisConfiguration* conf,
                    AnalysisConfiguration* ac,
                    EventFilter* eventFilter,
                    Event* event,
                    Event* recoEvent,
                    Event* uncorrEvent)
{
  /* for having the balance function correctly extracted the particle filters have to follow certain order */
  /* - charged particle should come always first                                                           */
  /* - particles of the same species, the positive has to come first and the negative immediately after    */
  /* - after the charged particles the might come any number of neutral                                    */
  /* - the balance function produced for neutrals will not have any sense                                  */
  std::vector<ParticleFilter<r>*> particleFilters;
  for (auto& part : conf->tpairs) {
    auto filter = PythiaAnalysisConfiguration::particleFilter<r>(part, efd, ac);
    if (filter != nullptr) {
      particleFilters.push_back(filter);
    } else {
      return false;
    }
  }
  /* the pairs taskname */
  TString taskName = TString::Format(conf->taskname.c_str(), TString::Format("PairsFDRej%s", efd.c_str()).Data());

  /* the two-particle mixed-event analyzer (raw pass) */
  analysisTasks[iTask++] = new TwoPartDiffCorrelationAnalyzerME<r, options>(taskName, ac, event, eventFilter, particleFilters);

  /* single particle analysis filters and task if any (raw pass) */
  if (conf->tsingles.size() > 0) {
    int nParticleFilters = 0;
    TString singlesTtaskName = TString::Format(conf->taskname.c_str(), TString::Format("SinglesFDRej%s", efd.c_str()).Data());
    ParticleFilter<r>** singleParticleFilters = new ParticleFilter<r>*[50];
    for (auto& part : conf->tsingles) {
      auto filter = PythiaAnalysisConfiguration::particleFilter<r>(part, efd, ac);
      if (filter != nullptr) {
        singleParticleFilters[nParticleFilters++] = filter;
      } else {
        return false;
      }
    }
    analysisTasks[iTask++] = new ParticleAnalyzer<r>(singlesTtaskName, ac, event, eventFilter, nParticleFilters, singleParticleFilters);
  }

  /* detector-effects passes: the corrected ("DetCorr", weights 1/ε, reading     */
  /* recoEvent) and the uncorrected ("Det", weights 1.0, reading uncorrEvent)    */
  /* analyzers. Both reuse the raw analyzer classes with their own fresh filter  */
  /* objects; the "DetCorr"/"Det" infix in the task name makes                   */
  /* Task::saveHistograms put each pass in its own TDirectory in the shared      */
  /* output .root file, side-by-side with the raw directories.                   */
  auto addDetectorPass = [&](const char* infix, Event* passEvent) -> bool {
    std::vector<ParticleFilter<r>*> passFilters;
    for (auto& part : conf->tpairs) {
      auto filter = PythiaAnalysisConfiguration::particleFilter<r>(part, efd, ac);
      if (filter == nullptr) {
        return false;
      }
      passFilters.push_back(filter);
    }
    TString passTaskName = TString::Format(conf->taskname.c_str(), TString::Format("Pairs%sFDRej%s", infix, efd.c_str()).Data());
    analysisTasks[iTask++] = new TwoPartDiffCorrelationAnalyzerME<r, options>(passTaskName, ac, passEvent, eventFilter, passFilters);

    if (conf->tsingles.size() > 0) {
      int nPassSingleFilters = 0;
      TString passSinglesTaskName = TString::Format(conf->taskname.c_str(), TString::Format("Singles%sFDRej%s", infix, efd.c_str()).Data());
      ParticleFilter<r>** passSingleFilters = new ParticleFilter<r>*[50];
      for (auto& part : conf->tsingles) {
        auto filter = PythiaAnalysisConfiguration::particleFilter<r>(part, efd, ac);
        if (filter == nullptr) {
          return false;
        }
        passSingleFilters[nPassSingleFilters++] = filter;
      }
      analysisTasks[iTask++] = new ParticleAnalyzer<r>(passSinglesTaskName, ac, passEvent, eventFilter, nPassSingleFilters, passSingleFilters);
    }
    return true;
  };

  if (conf->detectoreffects && recoEvent != nullptr && uncorrEvent != nullptr) {
    if (!addDetectorPass("DetCorr", recoEvent)) {
      return false;
    }
    if (!addDetectorPass("Det", uncorrEvent)) {
      return false;
    }
  }
  return true;
}

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
  std::vector<float> ptRangeLows = conf->ptRangeLows;
  std::vector<float> ptRangeUps = conf->ptRangeUps;
  std::vector<int> nPtRangeBins = conf->nPtRangeBins;
  if (ptRangeLows.size() != ptRangeUps.size() || ptRangeLows.size() != nPtRangeBins.size()) {
    Fatal("RunPythiaSimulationTwoParticlesDiff", "The limits or the number of bins of the pT ranges differ. Please, fix it!");
  }

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
  /* parallel events for the detector-effects passes: recoEvent holds the        */
  /* corrected ("DetCorr") particle set (weights 1/ε), uncorrEvent holds the     */
  /* same set uncorrected ("Det", weights 1.0). Both created iff detectoreffects.*/
  Event* recoEvent = nullptr;
  Event* uncorrEvent = nullptr;
  /* effHistos is parallel to conf->tpairs: effHistos[i] is the efficiency TH1 */
  /* for particles tagged with ixID == i (i.e. accepted by tpairFilters[i]).   */
  /* Missing entries (nullptr) => eff=1 for that index (silent fallback).      */
  std::vector<TH1*> effHistos;
  if (!conf->inputfile.empty()) {
    TFile* f = new TFile(conf->inputfile.c_str());
    if (f != nullptr && f->IsOpen()) {
      bool status = TH1::AddDirectoryStatus();
      TH1::AddDirectory(false);
      event->setMultiplicityPercentiles(f);
      /* load `<tpairs[i]>Efficiency` for each tpairs entry; missing -> nullptr */
      if (conf->detectoreffects) {
        effHistos.assign(conf->tpairs.size(), nullptr);
        for (size_t i = 0; i < conf->tpairs.size(); ++i) {
          std::string hname = conf->tpairs[i] + "Efficiency";
          TObject* o = f->Get(hname.c_str());
          if (o == nullptr) {
            Warning("RunPythiaSimulationTwoParticlesDiffME",
                    "efficiency histogram '%s' not found in %s -- eff=1 for tpairs[%zu]=%s",
                    hname.c_str(), conf->inputfile.c_str(), i, conf->tpairs[i].c_str());
            continue;
          }
          effHistos[i] = (TH1*)o->Clone(); /* detached: survives f->Close() */
        }
      }
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
  /* create the parallel reconstructed Events used by the detector-effects passes */
  if (conf->detectoreffects) {
    if (effHistos.empty()) {
      /* no input file or no histos loaded; keep effHistos sized to tpairs (all nullptr) */
      effHistos.assign(conf->tpairs.size(), nullptr);
    }
    recoEvent = Event::createReconstructed();
    uncorrEvent = Event::createReconstructed();
  }

  // ==========================
  // Generator Section
  // ==========================
  float genMinPt = conf->genMinPt;
  float genMaxPt = conf->genMaxPt;
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
  Task* detTask = nullptr;
  /* throwaway minimal configuration for the DetectorEffectsTask: all lifecycle flags */
  /* are false so Task::initialize / reset / finalize do nothing (the task produces   */
  /* no output of its own; the reconstructed Event is consumed by parallel analyzers).*/
  AnalysisConfiguration* detTaskCfg = nullptr;
  if (conf->detectoreffects) {
    detTaskCfg = new AnalysisConfiguration("DETEFF", "DETEFF", "1.0");
    detTaskCfg->loadHistograms = false;
    detTaskCfg->createHistograms = false;
    detTaskCfg->scaleHistograms = false;
    detTaskCfg->calculateDerivedHistograms = false;
    detTaskCfg->saveHistograms = false;
    detTaskCfg->resetHistograms = false;
    detTaskCfg->clearHistograms = false;
    detTaskCfg->forceHistogramsRewrite = false;
  }
  /* particle selection at the generator level */
  if (conf->gparticlefilter == "All" && conf->gchargefilter == "All") {
    if (conf->inrapidity) {
      ParticleFilter<AnalysisConfiguration::kRapidity>* particleFilterGen = new ParticleFilter<AnalysisConfiguration::kRapidity>(ParticleFilter<AnalysisConfiguration::kRapidity>::AllSpecies,
                                                                                                                                 ParticleFilter<AnalysisConfiguration::kRapidity>::AllCharges,
                                                                                                                                 ParticleFilter<AnalysisConfiguration::kRapidity>::None,
                                                                                                                                 genMinPt, genMaxPt,
                                                                                                                                 -abs_y[0], abs_y[0]);
      generator = new PythiaEventGenerator<AnalysisConfiguration::kRapidity>("PYTHIA", pc, event, eventFilterGen, particleFilterGen);
      if (conf->detectoreffects) {
        detTask = new DetectorEffectsTask<AnalysisConfiguration::kRapidity>("DETEFFECTS", detTaskCfg, event, recoEvent, uncorrEvent,
                                                                            conf->tpairs, effHistos,
                                                                            conf->mergedeta, conf->mergedphi, conf->mergedpt, seed);
      }
    } else {
      ParticleFilter<AnalysisConfiguration::kPseudorapidity>* particleFilterGen = new ParticleFilter<AnalysisConfiguration::kPseudorapidity>(ParticleFilter<AnalysisConfiguration::kPseudorapidity>::AllSpecies,
                                                                                                                                             ParticleFilter<AnalysisConfiguration::kPseudorapidity>::AllCharges,
                                                                                                                                             ParticleFilter<AnalysisConfiguration::kPseudorapidity>::None,
                                                                                                                                             genMinPt, genMaxPt,
                                                                                                                                             -abs_y[0], abs_y[0]);
      generator = new PythiaEventGenerator<AnalysisConfiguration::kPseudorapidity>("PYTHIA", pc, event, eventFilterGen, particleFilterGen);
      if (conf->detectoreffects) {
        detTask = new DetectorEffectsTask<AnalysisConfiguration::kPseudorapidity>("DETEFFECTS", detTaskCfg, event, recoEvent, uncorrEvent,
                                                                                  conf->tpairs, effHistos,
                                                                                  conf->mergedeta, conf->mergedphi, conf->mergedpt, seed);
      }
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

  analysisTasks = new Task*[nAnalysisTasks];
  iTask = 0;

  for (float y : abs_y) {
    double min_y = -y;
    double max_y = y;
    int nBins_y = int((max_y - min_y) / 0.1);

    for (uint iPtRange = 0; iPtRange < ptRangeLows.size(); ++iPtRange) {

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
      ac->rootOuputFileName = TString::Format("%s_%03d", TString::Format(conf->outputfname.c_str(), int(y * 10), int(ptRangeLows[iPtRange] * 10), int(ptRangeUps[iPtRange] * 10)).Data(), jobix).Data();
      ac->outputDirectory = TString::Format("%s", TString::Format(conf->outputfname.c_str(), int(y * 10), int(ptRangeLows[iPtRange] * 10), int(ptRangeUps[iPtRange] * 10)).Data()).Data();
      ac->histoBaseName = "TEST";

      ac->bin_edges_pt = conf->ptbins;
      ac->nBins_pt = nPtRangeBins[iPtRange];
      ac->min_pt = ptRangeLows[iPtRange];
      ac->max_pt = ptRangeUps[iPtRange];
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
      ac->fillPratt = conf->fillpratt;
      ac->fillInvariantMass = conf->fillinvmass;
      if (conf->inrapidity) {
        ac->fillYorEta = AnalysisConfiguration::kRapidity;
      } else {
        ac->fillYorEta = AnalysisConfiguration::kPseudorapidity;
      }

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
        for (auto& fd : conf->tfeeddownrej) {
          /* for having the balance function correctly extracted the particle filters have to follow certain order */
          /* - charged particle should come always first                                                           */
          /* - particles of the same species, the positive has to come first and the negative immediately after    */
          /* - after the charged particles the might come any number of neutral                                    */
          /* - the balance function produced for neutrals will not have any sense                                  */
          if (conf->inrapidity) {
            if (conf->fillpratt || conf->fillinvmass) {
              if (conf->fillinvmass) {
                if (conf->fillpratt) {
                  if (!configureTasks<AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillPrattAndInvariantMass>(fd, conf, ac, eventFilter, event, recoEvent, uncorrEvent)) {
                    return 0;
                  }
                } else {
                  if (!configureTasks<AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillInvariantMass>(fd, conf, ac, eventFilter, event, recoEvent, uncorrEvent)) {
                    return 0;
                  }
                }
              } else {
                if (!configureTasks<AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillPratt>(fd, conf, ac, eventFilter, event, recoEvent, uncorrEvent)) {
                  return 0;
                }
              }
            } else {
              if (!configureTasks<AnalysisConfiguration::kRapidity, AnalysisConfiguration::kNoAdditionalOptions>(fd, conf, ac, eventFilter, event, recoEvent, uncorrEvent)) {
                return 0;
              }
            }
          } else {
            if (conf->fillpratt || conf->fillinvmass) {
              if (conf->fillinvmass) {
                if (conf->fillpratt) {
                  if (!configureTasks<AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillPrattAndInvariantMass>(fd, conf, ac, eventFilter, event, recoEvent, uncorrEvent)) {
                    return 0;
                  }
                } else {
                  if (!configureTasks<AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillInvariantMass>(fd, conf, ac, eventFilter, event, recoEvent, uncorrEvent)) {
                    return 0;
                  }
                }
              } else {
                if (!configureTasks<AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillPratt>(fd, conf, ac, eventFilter, event, recoEvent, uncorrEvent)) {
                  return 0;
                }
              }
            } else {
              if (!configureTasks<AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kNoAdditionalOptions>(fd, conf, ac, eventFilter, event, recoEvent, uncorrEvent)) {
                return 0;
              }
            }
          }
        }
      }
    }
  }

  nAnalysisTasks = iTask;

  // ==========================
  // Event Loop
  // ==========================

  EventLoop* eventLoop = new EventLoop();
  generator->reportLevel = repLevel;
  eventLoop->addTask((Task*)generator);
  /* the DetectorEffectsTask (if any) MUST run after the generator and before the */
  /* reconstructed analyzers; it builds the parallel recoEvent each event.        */
  if (detTask != nullptr) {
    detTask->reportLevel = repLevel;
    eventLoop->addTask(detTask);
  }
  for (int iAnalysisTask = 0; iAnalysisTask < nAnalysisTasks; iAnalysisTask++) {
    analysisTasks[iAnalysisTask]->reportLevel = repLevel;
    eventLoop->addTask(analysisTasks[iAnalysisTask]);
  }
  eventLoop->run(nEventsRequested, nEventsReport);

  TH1::AddDirectory(oldstatus);
  cout << "<INFO> PYTHIA Model Analysis - Pair Differential Correlations Histograms  - Completed" << endl;
  time(&end); // note time after execution
  double difference = difftime(end, begin);
  cout << "<INFO> in " << difference << " seconds" << endl;
}
