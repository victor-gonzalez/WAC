#ifndef PYTHIAANALYSISCONFIGURATION_H
#define PYTHIAANALYSISCONFIGURATION_H

#include <TObject.h>
#include <TError.h>
#include "ParticleFilter.hpp"
#include "AnalysisConfiguration.hpp"

class PythiaAnalysisConfiguration : public TObject
{
 public:
  PythiaAnalysisConfiguration() {}
  ~PythiaAnalysisConfiguration() {}

 public:
  std::vector<float> abs_y = {10.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.8};
  float min_pt = 0.0;
  float max_pt = 1e6;
  std::vector<double> ptbins = {0.0,
                                0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                                0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4,
                                1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                                3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 13.0, 20.0};
  int n_ptbins = ptbins.size() - 1;

  long nEventsRequested = 500000;
  int nEventsReport = 10000;
  std::string logLevel = "info";
  std::vector<std::string> pythiaOptions = {
    "Init:showChangedSettings = on",
    "Init:showChangedParticleData = on",
    "Next:numberCount = 10000",
    "Next:numberShowInfo = 1",
    "Next:numberShowProcess = 0",
    "Next:numberShowEvent = 0",
    "SoftQCD:all = on",
    "Tune:pp = 14"};
  int projectileA = 2212;
  int projectileB = 2212;
  float energy = 7000.0;
  std::string geventfilter = "MB";
  std::string gparticlefilter = "All";
  std::string gchargefilter = "All";
  bool inrapidity = true;
  std::string outputfname = "PYTHIA8_Rapidity%03d";
  std::string taskname = "%s_All";
  std::string inputfile = "MultPercentile.root";
  std::vector<std::string> teventfilter = {"MB"};
  std::vector<std::string> tpairs = {"AllP", "AllM", "La", "ALa", "Gam"};
  std::vector<std::string> tsingles = {"AllA"};
  std::vector<std::string> tfeeddownrej = {"none", "all"};

  template <AnalysisConfiguration::RapidityPseudoRapidity r>
  static ParticleFilter<r>* particleFilter(std::string str, std::string fdstr, AnalysisConfiguration* ac)
  {
    auto getparticle = [](auto str) {
      if (str == "PiP" || str == "PiM" || str == "PiC" || str == "Pi0" || str == "PiA") {
        return ParticleFilter<r>::Pion;
      } else if (str == "KaP" || str == "KaM" || str == "KaC" || str == "Ka0" || str == "KaA") {
        return ParticleFilter<r>::Kaon;
      } else if (str == "PrP" || str == "PrM" || str == "PrC" || str == "PrA") {
        return ParticleFilter<r>::Proton;
      } else if (str == "La") {
        return ParticleFilter<r>::Lambda;
      } else if (str == "ALa") {
        return ParticleFilter<r>::ALambda;
      } else if (str == "Gam") {
        return ParticleFilter<r>::Photon;
      } else if (str == "AllP" || str == "AllM" || str == "AllC" || str == "All0" || str == "AllA") {
        return ParticleFilter<r>::AllSpecies;
      } else {
        ::Fatal("PythiaAnalysisConfiguration::particleFilter()", "Paricle species %s still not supported for analysis. Please fix it!!", str.c_str());
        return ParticleFilter<r>::AllSpecies;
      }
    };
    auto getcharge = [](auto str) {
      if (str == "PiP" || str == "KaP" || str == "PrP" || str == "AllP") {
        return ParticleFilter<r>::Positive;
      } else if (str == "PiM" || str == "KaM" || str == "PrM" || str == "AllM") {
        return ParticleFilter<r>::Negative;
      } else if (str == "PiC" || str == "KaC" || str == "PrC" || str == "AllC") {
        return ParticleFilter<r>::Charged;
      } else if (str == "Pi0" || str == "Ka0" || str == "All0" || str == "La" || str == "ALa" || str == "Gam") {
        return ParticleFilter<r>::Neutral;
      } else if (str == "PiA" || str == "KaA" || str == "PrA" || str == "AllA") {
        return ParticleFilter<r>::AllCharges;
      } else {
        ::Fatal("PythiaAnalysisConfiguration::particleFilter()", "Paricle species %s still not supported for analysis. Please fix it!!", str.c_str());
        return ParticleFilter<r>::AllCharges;
      }
    };
    auto getfeeddown = [](auto str) {
      if (str == "none") {
        return ParticleFilter<r>::None;
      } else if (str == "all") {
        return ParticleFilter<r>::AllResonances;
      } else {
        ::Fatal("PythiaAnalysisConfiguration::particleFilter()", "Resonances suppression %s still not supported for analysis. Please fix it!!", str.c_str());
        return ParticleFilter<r>::None;
      }
    };
    return new ParticleFilter<r>(getparticle(str), getcharge(str), getfeeddown(fdstr), ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
  }

  ClassDef(PythiaAnalysisConfiguration, 2)
};

#endif // PYTHIAANALYSISCONFIGURATION_H
