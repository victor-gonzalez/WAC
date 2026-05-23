// Author: Victor Gonzalez, 2026

/***********************************************************************
 * Copyright (C) 2026.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#include <vector>
#include "TFile.h"
#include "TMath.h"
#include "DetectorEffectsTask.hpp"
#include "Particle.hpp"

templateClassImp(DetectorEffectsTask);

template <AnalysisConfiguration::RapidityPseudoRapidity r>
DetectorEffectsTask<r>::DetectorEffectsTask(const TString& name,
                                            TaskConfiguration* configuration,
                                            Event* _event,
                                            Event* _dstEvent,
                                            Event* _dstUncorrEvent,
                                            const std::vector<std::string>& _trackNames,
                                            const std::vector<TH1*>& _effHistos,
                                            double _dEta,
                                            double _dPhi,
                                            double _dPt,
                                            long rngSeed,
                                            MergeHandling _mergeMode)
  : Task(name, configuration, _event),
    dstEvent(_dstEvent),
    dstUncorrEvent(_dstUncorrEvent),
    trackNames(_trackNames),
    effHistos(_effHistos),
    dEta(_dEta),
    dPhi(_dPhi),
    dPt(_dPt),
    mergeMode(_mergeMode),
    rng(new TRandom3(rngSeed)),
    hEvtCount_vs_nSrc(nullptr),
    hNAfterEffSum_vs_nSrc(nullptr),
    hNDropMergeSum_vs_nSrc(nullptr),
    profFracMerged_vs_nSrc(nullptr)
{
  if (!event) {
    if (reportError())
      cout << "DetectorEffectsTask::CTOR(...) src Event is a null pointer." << endl;
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
  /* the QA histograms are not registered with any TDirectory (we call    */
  /* SetDirectory(0) at allocation), so they are owned exclusively by us. */
  delete hEvtCount_vs_nSrc;
  delete hNAfterEffSum_vs_nSrc;
  delete hNDropMergeSum_vs_nSrc;
  delete profFracMerged_vs_nSrc;
}

/// \brief Allocate the QA histograms
template <AnalysisConfiguration::RapidityPseudoRapidity r>
void DetectorEffectsTask<r>::createHistograms()
{
  const int    nbins = 300;
  const double xlo   = 0.0;
  const double xhi   = 3000.0; /* generous upper edge -- Pythia pp events sit well below */

  const TString base = getName();
  hEvtCount_vs_nSrc = new TH1F(base + "_hEvtCount_vs_nSrc",
                               "events vs n_{src};n_{src};events",
                               nbins, xlo, xhi);
  hEvtCount_vs_nSrc->SetDirectory(0);
  hNAfterEffSum_vs_nSrc = new TProfile(base + "_hNAfterEffSum_vs_nSrc",
                                       "#Sigma n_{afterEff} vs n_{src};n_{src};#Sigma n_{afterEff}",
                                       nbins, xlo, xhi);
  hNAfterEffSum_vs_nSrc->SetDirectory(0);
  hNDropMergeSum_vs_nSrc = new TProfile(base + "_hNDropMergeSum_vs_nSrc",
                                        "#Sigma n_{dropped by merging} vs n_{src};n_{src};#Sigma n_{drop}",
                                        nbins, xlo, xhi);
  hNDropMergeSum_vs_nSrc->SetDirectory(0);
  profFracMerged_vs_nSrc = new TProfile(base + "_profFracMerged_vs_nSrc",
                                        "fraction merged vs n_{src};n_{src};#LT n_{drop}/n_{afterEff} #GT",
                                        nbins, xlo, xhi);
  profFracMerged_vs_nSrc->SetDirectory(0);
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
  if (!event) {
    if (reportError())
      cout << "DetectorEffectsTask::execute() src event is null. Abort." << endl;
    postTaskError();
    return;
  }

  if (dstEvent == nullptr) {
    /* nothing to do */
    return;
  }

  /* clear the reconstructed event and inherit its centrality / multiplicity */
  /* class from the raw event (so EventFilter::accept and the analyzers'     */
  /* fillEventWiseInfo see consistent values)                                */
  dstEvent->reset();
  dstEvent->copyEventLevelInfoFrom(*event);

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
  const int nSrc = int(event->getNParticles());
  int nKept = 0;
  for (int i = 0; i < nSrc; ++i) {
    Particle* src = event->getParticleAt(i);
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
  const int nAfterEff = nKept; /* captured BEFORE the merging stage for QA */

  /* ============================================================ */
  /* Stage 2: track merging  (currently only kKeepLeading wired)  */
  /*                                                              */
  /* SAME-CHARGE gate: in a magnetic-field tracker opposite-charge */
  /* tracks curve apart and stop sharing hits quickly, so the      */
  /* merging channel is dominantly same-sign; opposite-sign close  */
  /* pairs are skipped here.                                       */
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
        if (a.charge != b.charge)
          continue; /* opposite-charge close pairs do not merge in real trackers */
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
  dstUncorrEvent->copyEventLevelInfoFrom(*event);
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

  /* ============================================================ */
  /* QA fills (always on)                                          */
  /* ============================================================ */
  const int nDropMerge = nAfterEff - nKept;
  if (hEvtCount_vs_nSrc != nullptr) {
    hEvtCount_vs_nSrc->Fill(nSrc);
    hNAfterEffSum_vs_nSrc->Fill(nSrc, nAfterEff);
    hNDropMergeSum_vs_nSrc->Fill(nSrc, nDropMerge);
    if (nAfterEff > 0)
      profFracMerged_vs_nSrc->Fill(nSrc, double(nDropMerge) / double(nAfterEff));
  }

  if (reportDebug())
    cout << "DetectorEffectsTask::execute() nGen=" << nSrc
         << " nAfterEff=" << nAfterEff
         << " nReco=" << nKept << endl;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void DetectorEffectsTask<r>::saveHistograms(TDirectory* dir)
{
  if (reportDebug())
    cout << "DetectorEffectsTask::saveHistograms(...) Saving Event histograms to file." << endl;
  if (!dir) {
    if (reportError())
      cout << "DetectorEffectsTask::saveHistograms(...) output directory is a null  pointer." << endl;
    postTaskError();
    return;
  }
  dir->cd();

  /* now save the event histograms */
  if (reportDebug())
    cout << "TwoPartDiffCorrelationAnalyzer::saveHistograms(...) saving event histograms." << endl;
  event->saveHistograms(dir);

  /* now save the detector effect histograms */
  if (reportDebug())
    cout << "DetectorEffectsTask::saveHistograms(...) saving event histograms." << endl;

  if (hEvtCount_vs_nSrc == nullptr) {
    if (reportError())
      cout << "DetectorEffectsTask::finalize() QA histograms not allocated; nothing to write." << endl;
    return;
  }

  dir->cd();
  hEvtCount_vs_nSrc->Write();
  hNAfterEffSum_vs_nSrc->Write();
  hNDropMergeSum_vs_nSrc->Write();
  profFracMerged_vs_nSrc->Write();
}

template class DetectorEffectsTask<AnalysisConfiguration::kRapidity>;
template class DetectorEffectsTask<AnalysisConfiguration::kPseudorapidity>;

templateClassImp(DetectorEffectsTask)
