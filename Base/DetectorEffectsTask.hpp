// Author: Victor Gonzalez, 2026

/***********************************************************************
 * Copyright (C) 2026.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_DetectorEffectsTask
#define WAC_DetectorEffectsTask

#include <string>
#include <vector>
#include "TH1.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "Task.hpp"
#include "Event.hpp"
#include "ParticleFilter.hpp"
#include "AnalysisConfiguration.hpp"

///////////////////////////////////////////////////////////////////////////////
// DetectorEffectsTask<r>
//
// Task inserted in the EventLoop between the event generator and the analysis
// tasks of the detector-effects passes.  Reads the singleton (raw) Event
// filled by the generator and produces TWO parallel Event instances, each
// carrying its own Factory<Particle> (created via Event::createReconstructed()):
//
//   * the CORRECTED ("DetCorr") event -- detector effects applied, every
//     surviving particle carrying the efficiency-correction weight 1/ε;
//   * the UNCORRECTED ("Det") event -- the very same reconstructed particle
//     set, but with every weight reset to 1.0 (raw detector response).
//
// The corrected event is built first by:
//
//   1) Efficiency stage. Each raw particle is classified by running it
//      through ParticleFilter<r>::getIndex(trackNames, p), which iterates
//      the `tpairs` track-name strings (handed in from the analysis
//      configuration) and returns the first identity-matching index
//      (species + charge only; kinematic and feedDown cuts are ignored,
//      so the classification is independent of the per-analyzer
//      acceptance windows that differ across the EventLoop's analyzers).
//      The efficiency ε(pT) for that index is read from `effHistos[ixID]`
//      (parallel to tpairs); the particle is dropped with probability
//      (1 − ε), and survivors are copied to the reconstructed Event with
//      weight = 1/ε.  Particles matching no track name, or whose
//      `effHistos[ixID]` is null, survive with ε = 1 and weight = 1.0.
//
//   2) Track-merging stage. For every pair of SAME-CHARGE surviving
//      particles within the configured proximity window
//      (|Δη|, |Δφ|, |ΔpT|), apply the merge handling: kKeepLeading drops
//      the lower-pT track of each too-close pair (keeping the higher-pT
//      one unchanged).  Same-sign gating is physical: in a magnetic-field
//      tracker, opposite-charge tracks curve apart and stop sharing hits
//      quickly, so merging is dominantly same-sign.  Other modes are
//      stubbed for a future change; only kKeepLeading is wired right now.
//      Merging removes both pairs and singles from the reconstructed
//      set; no correction is applied for it.
//
//   3) Uncorrected replication.  Once the corrected event is complete it is
//      copied particle-by-particle into the uncorrected event with every
//      weight reset to 1.0, so the parallel "Det" analyzers see the detector
//      response without the efficiency correction.
//
// Per-event RNG: the task owns its own TRandom3 seeded from the constructor
// argument (the same per-job `seed` used by Pythia), so the new draws are
// independent across jobs without disturbing gRandom (which the raw pass
// relies on for the event-plane angle).
//
// Templated on r (kRapidity / kPseudorapidity) for symmetry with the
// analyzers.  The proximity test uses pseudo-rapidity (η) for both builds:
// the detector geometry is η-based, regardless of which variable the
// analysis bins in.
//
// QA histograms (always-on, self-contained -- the task ignores the lifecycle
// flags on its TaskConfiguration): per-event counts binned in n_src (number
// of generated particles in the raw Event) are accumulated into a small set
// of TH1 / TProfile and written, when saveHistograms() runs.
///////////////////////////////////////////////////////////////////////////////

template <AnalysisConfiguration::RapidityPseudoRapidity r>
class DetectorEffectsTask : public Task
{
 public:
  /// How a too-close pair should be handled.
  enum MergeHandling {
    kKeepLeading = 0,    ///< drop the lower-pT track, keep the higher-pT one unchanged (CURRENT BEHAVIOUR)
    kRemoveBoth = 1,     ///< drop both tracks (stub; not wired to a config switch yet)
    kCollapseToOne = 2   ///< replace the two by one combined track (stub; not wired)
  };

  DetectorEffectsTask(const TString& name,
                      TaskConfiguration* configuration,
                      Event* srcEvent,
                      Event* dstRecoEvent,
                      Event* dstUncorrEvent,
                      const std::vector<std::string>& trackNames,
                      const std::vector<TH1*>& effHistos,
                      double dEta,
                      double dPhi,
                      double dPt,
                      long rngSeed,
                      MergeHandling mergeMode = kKeepLeading);
  virtual ~DetectorEffectsTask();

  virtual void createHistograms();
  virtual void execute();
  virtual void saveHistograms(TDirectory* dir);

 private:
  double efficiencyForIxID(int ixID, double pt) const;
  static double wrapDeltaPhi(double dphi);
  bool tooClose(double etaA, double phiA, double ptA,
                double etaB, double phiB, double ptB) const;
  bool mergeActive() const { return dEta > 0.0 && dPhi > 0.0 && dPt > 0.0; }

  Event* dstEvent;                     ///< corrected ("DetCorr") event; surviving particles carry weight 1/ε
  Event* dstUncorrEvent;               ///< uncorrected ("Det") event; same particle set with weight 1.0
  std::vector<std::string> trackNames; ///< the `tpairs` vocabulary, parallel to effHistos
  std::vector<TH1*> effHistos;         ///< parallel to trackNames; nullptr entries => ε=1
  double dEta;
  double dPhi;
  double dPt;
  MergeHandling mergeMode;
  TRandom3* rng; ///< task-owned RNG (does not touch gRandom)

  /* QA: per-event distributions binned in n_src (number of raw particles). */
  /* hEvtCount_vs_nSrc is the event-count denominator; the two _Sum_ TH1s   */
  /* hold sum-of-nAfterEff and sum-of-nDropped(merging) per nSrc bin -- the */
  /* ratio of the two yields an hadd-unbiased fraction-merged-per-event.    */
  /* profFracMerged_vs_nSrc is the same information in TProfile form, more  */
  /* convenient to look at directly.                                        */
  TH1F* hEvtCount_vs_nSrc;
  TProfile* hNAfterEffSum_vs_nSrc;
  TProfile* hNDropMergeSum_vs_nSrc;
  TProfile* profFracMerged_vs_nSrc;

  ClassDef(DetectorEffectsTask, 0)
};

#endif /* WAC_DetectorEffectsTask */
