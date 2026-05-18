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
  std::vector<float> ptRangeLows = {0.2, 0.2, 0.2, 0.2};
  std::vector<float> ptRangeUps = {0.6, 1.0, 1.5, 2.0};
  std::vector<int> nPtRangeBins = {4, 8, 13, 18};
  std::vector<double> ptbins = {0.0,
                                0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                                0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4,
                                1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                                3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 13.0, 20.0};

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
  float genMinPt = 0.2;
  float genMaxPt = 2.0;
  std::string geventfilter = "MB";
  std::string gparticlefilter = "All";
  std::string gchargefilter = "All";
  bool inrapidity = true;
  bool fillpratt = false;
  bool fillinvmass = false;
  std::string outputfname = "PYTHIA8_Rapidity%03d_Pt%02d%02d";
  std::string taskname = "%s_All";
  std::string inputfile = "MultPercentile.root";
  std::vector<std::string> teventfilter = {"MB"};
  std::vector<std::string> tpairs = {"AllP", "AllM", "La", "ALa", "Gam"};
  std::vector<std::string> tsingles = {"AllA"};
  std::vector<std::string> tfeeddownrej = {"none", "all"};

  /* ============================================================== */
  /* Detector effects (second analysis pass).                       */
  /* `detectoreffects` is the master enable: when true, the         */
  /* launcher inserts a DetectorEffectsTask after the generator and */
  /* duplicates each analyzer onto a parallel reconstructed Event.  */
  /* Within the reco pass each effect is independently switchable:  */
  /*                                                                */
  /*  - Efficiency drop: keyed by the index of the particle in the  */
  /*    `tpairs` vector (the same ixID the analyzers assign).  For  */
  /*    each entry `pname` in `tpairs` the input ROOT file is       */
  /*    expected to contain a TH1 named `<pname>Efficiency` whose   */
  /*    x-axis is pT [GeV/c] and whose bin contents are eff in      */
  /*    [0,1].  Missing histograms => eff=1 for that index (no      */
  /*    drop), with a warning printed at startup.  Particles that   */
  /*    are not accepted by any `tpairs` filter also survive with   */
  /*    eff=1, weight 1.0.                                          */
  /*                                                                */
  /*  - Track merging: active iff all three thresholds are >0.      */
  /*    Two tracks within (|Deta|, |Dphi|, |DpT|) are merged by     */
  /*    dropping the lower-pT one (keep-leading).  No correction    */
  /*    is applied for this effect for now.                         */
  /*                                                                */
  /* `detectoreffects:true` with no efficiency histograms in the    */
  /* input file and all thresholds <=0 => reco pass = exact copy of */
  /* raw (regression check).                                        */
  /* ============================================================== */
  bool detectoreffects = false;
  float mergedeta = -1.0;
  float mergedphi = -1.0;
  float mergedpt = -1.0;

  template <AnalysisConfiguration::RapidityPseudoRapidity r>
  static ParticleFilter<r>* particleFilter(std::string str, std::string fdstr, AnalysisConfiguration* ac)
  {
    /* (species, charge, feeddown) lookup tables now live on ParticleFilter; this   */
    /* method only adds the per-(y,pT) kinematic window from the supplied `ac`.     */
    return new ParticleFilter<r>(ParticleFilter<r>::speciesFor(str),
                                 ParticleFilter<r>::chargeFor(str),
                                 ParticleFilter<r>::feedDownFor(fdstr),
                                 ac->min_pt, ac->max_pt, ac->min_y, ac->max_y);
  }

  ClassDef(PythiaAnalysisConfiguration, 4)
};

#endif // PYTHIAANALYSISCONFIGURATION_H
