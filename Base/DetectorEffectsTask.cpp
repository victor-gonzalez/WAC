// Author: Victor Gonzalez, 2026

/***********************************************************************
 * Copyright (C) 2026.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#include <vector>
#include "TMath.h"
#include "DetectorEffectsTask.hpp"
#include "Particle.hpp"

templateClassImp(DetectorEffectsTask);

template <AnalysisConfiguration::RapidityPseudoRapidity r>
DetectorEffectsTask<r>::DetectorEffectsTask(const TString& name,
                                            TaskConfiguration* configuration,
                                            Event* _srcEvent,
                                            Event* _dstEvent,
                                            Event* _dstUncorrEvent,
                                            const std::vector<std::string>& _trackNames,
                                            const std::vector<TH1*>& _effHistos,
                                            double _dEta,
                                            double _dPhi,
                                            double _dPt,
                                            long rngSeed,
                                            MergeHandling _mergeMode)
  : Task(name, configuration, _srcEvent),
    srcEvent(_srcEvent),
    dstEvent(_dstEvent),
    dstUncorrEvent(_dstUncorrEvent),
    trackNames(_trackNames),
    effHistos(_effHistos),
    dEta(_dEta),
    dPhi(_dPhi),
    dPt(_dPt),
    mergeMode(_mergeMode),
    rng(new TRandom3(rngSeed))
{
  if (!srcEvent || !dstEvent || !dstUncorrEvent) {
    if (reportError())
      cout << "DetectorEffectsTask::CTOR(...) src or dst Event is a null pointer." << endl;
    postTaskError();
  }
  if (effHistos.size() != trackNames.size()) {
    if (reportError())
      cout << "DetectorEffectsTask::CTOR(...) effHistos size (" << effHistos.size()
           << ") differs from trackNames size (" << trackNames.size() << ")." << endl;
    postTaskError();
  }
  if (reportInfo()) {
    int nEffSet = 0;
    for (auto* h : effHistos)
      if (h != nullptr)
        ++nEffSet;
    cout << "DetectorEffectsTask::CTOR(...) " << getName()
         << " track-name entries: " << trackNames.size()
         << " of which with efficiency histo: " << nEffSet
         << " merging (dEta, dPhi, dPt) = (" << dEta << ", " << dPhi << ", " << dPt << ")"
         << (mergeActive() ? " [ACTIVE]" : " [inactive]") << endl;
  }
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
DetectorEffectsTask<r>::~DetectorEffectsTask()
{
  delete rng;
}

/// \brief Look up the detection efficiency for a particle.
///
/// Indexed by the `tpairs` ixID (the same internal ID the analyzers assign);
/// the launcher populates `effHistos` parallel to the `tpairs` vector by
/// loading the TH1 named `<tpairs[i]>Efficiency` from the input ROOT file.
/// If the ixID is negative (particle not accepted by any tpairs filter) or
/// `effHistos[ixID]` is null (histogram missing in the file), the
/// efficiency is 1.0 (the particle survives unchanged with weight 1.0).
template <AnalysisConfiguration::RapidityPseudoRapidity r>
double DetectorEffectsTask<r>::efficiencyForIxID(int ixID, double pt) const
{
  if (ixID < 0 || ixID >= int(effHistos.size()))
    return 1.0;
  const TH1* h = effHistos[ixID];
  if (h == nullptr)
    return 1.0;
  return h->GetBinContent(h->FindFixBin(pt));
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
double DetectorEffectsTask<r>::wrapDeltaPhi(double dphi)
{
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
bool DetectorEffectsTask<r>::tooClose(double etaA, double phiA, double ptA,
                                       double etaB, double phiB, double ptB) const
{
  /* detector geometry is in eta regardless of which variable the analysis bins in */
  if (TMath::Abs(etaA - etaB) >= dEta)
    return false;
  if (TMath::Abs(wrapDeltaPhi(phiA - phiB)) >= dPhi)
    return false;
  if (TMath::Abs(ptA - ptB) >= dPt)
    return false;
  return true;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void DetectorEffectsTask<r>::execute()
{
  if (!srcEvent || !dstEvent || !dstUncorrEvent) {
    if (reportError())
      cout << "DetectorEffectsTask::execute() src or dst Event is null. Abort." << endl;
    postTaskError();
    return;
  }

  /* clear the reconstructed event and inherit its centrality / multiplicity */
  /* class from the raw event (so EventFilter::accept and the analyzers'     */
  /* fillEventWiseInfo see consistent values)                                */
  dstEvent->reset();
  dstEvent->copyEventLevelInfoFrom(*srcEvent);

  /* ============================================================ */
  /* Stage 1: per-particle detection efficiency                   */
  /*                                                              */
  /* The classifying ixID is computed by running the particle     */
  /* through ParticleFilter<r>::getIndex(trackNames, *src) which  */
  /* matches on IDENTITY only (species + charge) -- kinematic     */
  /* and feedDown cuts are intentionally ignored so the           */
  /* classification is independent of any per-analyzer            */
  /* acceptance window.  Particles not matching any track name    */
  /* (ixID < 0) are kept unchanged with weight 1.0 so the         */
  /* merging stage can still see them.                            */
  /* ============================================================ */
  const int nSrc = int(srcEvent->getNParticles());
  int nKept = 0;
  for (int i = 0; i < nSrc; ++i) {
    Particle* src = srcEvent->getParticleAt(i);
    if (src == nullptr)
      continue;
    const int taggingIxID = trackNames.empty() ? -1
                                               : ParticleFilter<r>::getIndex(trackNames, *src);
    const double eps = efficiencyForIxID(taggingIxID, src->pt);
    if (eps < 1.0 && rng->Rndm() > eps)
      continue; /* dropped (eps==1.0 path skips the RNG draw for the bit-identical regression) */
    Particle* dst = dstEvent->appendParticle();
    if (dst == nullptr) {
      if (reportError())
        cout << "DetectorEffectsTask::execute() reconstructed event factory exhausted at nKept=" << nKept << endl;
      postTaskError();
      return;
    }
    *dst = *src;
    dst->weight = (eps > 0.0) ? 1.0 / eps : 0.0;
    dst->ixID = -1; /* the reco analyzers will recompute ixID via their own (narrower) filters */
    ++nKept;
  }

  /* expose the kept particles via getParticleAt so the merging stage can read them */
  dstEvent->setNParticlesAccepted(nKept);

  /* ============================================================ */
  /* Stage 2: track merging  (currently only kKeepLeading wired)  */
  /* ============================================================ */
  if (mergeActive() && nKept > 1) {
    std::vector<bool> dead(nKept, false);

    for (int i = 0; i < nKept; ++i) {
      if (dead[i])
        continue;
      Particle& a = *dstEvent->getParticleAt(i);
      for (int j = i + 1; j < nKept; ++j) {
        if (dead[j])
          continue;
        Particle& b = *dstEvent->getParticleAt(j);
        if (!tooClose(a.eta, a.phi, a.pt, b.eta, b.phi, b.pt))
          continue;
        /* keep-leading-drop-other: the lower-pT track of the pair dies */
        if (a.pt < b.pt) {
          dead[i] = true;
          break; /* i is dead now; move outer to the next i */
        } else {
          dead[j] = true;
          /* a may still be alive and may pair with another k > j */
        }
      }
    }

    /* compact survivors in place (w <= i always, so no aliasing hazard) */
    int w = 0;
    for (int i = 0; i < nKept; ++i) {
      if (dead[i])
        continue;
      if (w != i) {
        *dstEvent->getParticleAt(w) = *dstEvent->getParticleAt(i);
      }
      ++w;
    }
    dstEvent->setNParticlesAccepted(w);
    nKept = w;
  }

  /* ============================================================ */
  /* Stage 3: replicate the reconstructed set into the uncorrected */
  /* (Det) event -- same particles, detector effects applied, but  */
  /* with the efficiency-correction weight reset to 1.0 so the     */
  /* parallel "Det" analyzers see the raw detector response.       */
  /* ============================================================ */
  dstUncorrEvent->reset();
  dstUncorrEvent->copyEventLevelInfoFrom(*srcEvent);
  for (int i = 0; i < nKept; ++i) {
    Particle* corr = dstEvent->getParticleAt(i);
    Particle* unc = dstUncorrEvent->appendParticle();
    if (unc == nullptr) {
      if (reportError())
        cout << "DetectorEffectsTask::execute() uncorrected event factory exhausted at i=" << i << endl;
      postTaskError();
      return;
    }
    *unc = *corr;
    unc->weight = 1.0;
  }
  dstUncorrEvent->setNParticlesAccepted(nKept);

  if (reportDebug())
    cout << "DetectorEffectsTask::execute() nGen=" << nSrc << " nReco=" << nKept << endl;
}

template class DetectorEffectsTask<AnalysisConfiguration::kRapidity>;
template class DetectorEffectsTask<AnalysisConfiguration::kPseudorapidity>;

templateClassImp(DetectorEffectsTask)
