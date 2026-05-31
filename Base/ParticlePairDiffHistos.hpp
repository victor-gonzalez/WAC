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

 Class defining ParticlePairHistos
 */

#ifndef WAC_ParticlePairDiffHistos
#define WAC_ParticlePairDiffHistos

#include "Histograms.hpp"
#include "EventPool.hpp"

class ParticlePairDiffHistos : public Histograms
{
 public:
  ParticlePairDiffHistos(const TString& name,
                         AnalysisConfiguration* configuration,
                         LogLevel debugLevel);
  ParticlePairDiffHistos(TDirectory* dir,
                         const TString& name,
                         AnalysisConfiguration* configuration,
                         LogLevel debugLevel);
  virtual ~ParticlePairDiffHistos();
  void initialize();
  template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
  int getGlobalDeltaYEtaDeltaPhiBin(ParticleType1& p1, ParticleType2& p2);
  template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
  float getDeltaYEta(ParticleType1& particle1, ParticleType2& particle2);
  template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
  float getDeltaPhi(ParticleType1& particle1, ParticleType2& particle2);
  template <typename ParticleType1, typename ParticleType2>
  float getInvMass(ParticleType1& particle1, ParticleType2& particle2);
  template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
  void getPratt(ParticleType1& particle1, ParticleType2& particle2, double& kT, double& qinv, double& qlong, double& qside, double& qout);
  template <AnalysisConfiguration::RapidityPseudoRapidity r, AnalysisConfiguration::FillPairOptions options, typename ParticleType1, typename ParticleType2>
  void fill(ParticleType1& particle1, ParticleType2& particle2, double weight1, double weight2, double pTavg1 = 0.0, double pTavg2 = 0.0);
  template <AnalysisConfiguration::RapidityPseudoRapidity r, AnalysisConfiguration::FillPairOptions options, typename ParticleType1, typename ParticleType2>
  void fill(ParticleType1& particle1, ParticleType2& particle2, double weight1, double weight2, TH2* pTavg1, TH2* pTavg2);
  void loadHistograms(TDirectory* dir);
  void fillEventWiseInfo(float multiplicity, double nAccepted, float weight);

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  TH2* h_n2_ptPt;
  TProfile* p_n2_vsC;

  TH2* h_n2_DetaDphi;
  TProfile2D* p_n2_DetaDphi;
  TH2* h_ptpt_DetaDphi;
  TH2* h_dptdpt_DetaDphi;

  TH2* h_n2_DyDphi;
  TProfile2D* p_n2_DyDphi;
  TH2* h_ptpt_DyDphi;
  TH2* h_dptdpt_DyDphi;

  TH2* h_n2_QinvKt;
  TH2* h_n2_QlongKt;
  TH2* h_n2_QsideKt;
  TH2* h_n2_QoutKt;
  TH3* h_n2_QlongQsideQout;
  TH2* h_ptpt_QinvKt;
  TH2* h_ptpt_QlongKt;
  TH2* h_ptpt_QsideKt;
  TH2* h_ptpt_QoutKt;
  TH3* h_ptpt_QlongQsideQout;
  TH2* h_dptdpt_QinvKt;
  TH2* h_dptdpt_QlongKt;
  TH2* h_dptdpt_QsideKt;
  TH2* h_dptdpt_QoutKt;
  TH3* h_dptdpt_QlongQsideQout;

  TH1* h_invMass;
  TH2* h_Pi0GG_DetaDphi;
  TH2* h_Pi0GG_DyDphi;
  TH2* h_Pi0GGSide_DetaDphi;
  TH2* h_Pi0GGSide_DyDphi;

  TH2* h_EtaGG_DetaDphi;
  TH2* h_EtaGG_DyDphi;
  TH2* h_EtaGGSide_DetaDphi;
  TH2* h_EtaGGSide_DyDphi;

  TH2* h_PhiKK_DetaDphi;
  TH2* h_PhiKK_DyDphi;
  TH2* h_PhiKKSide_DetaDphi;
  TH2* h_PhiKKSide_DyDphi;

  ClassDef(ParticlePairDiffHistos, 6)
};

/// WARNING: for performance reasons no checks are done about the consistency
/// of particles' eta and phi within the corresponding ranges so, it is suppossed
/// the particles have been accepted and they are within that ranges
/// IF THAT IS NOT THE CASE THE ROUTINE WILL PRODUCE NONSENSE RESULTS
template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
inline float ParticlePairDiffHistos::getDeltaYEta(ParticleType1& particle1, ParticleType2& particle2)
{
  if constexpr (r == AnalysisConfiguration::kRapidity) {
    return h_n2_DyDphi->GetXaxis()->GetBinCenter(configuration->getDeltaYEtaIndex<r>(particle1, particle2) + 1);
  } else {
    return h_n2_DetaDphi->GetXaxis()->GetBinCenter(configuration->getDeltaYEtaIndex<r>(particle1, particle2) + 1);
  }
}

template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
inline float ParticlePairDiffHistos::getDeltaPhi(ParticleType1& particle1, ParticleType2& particle2)
{
  if constexpr (r == AnalysisConfiguration::kRapidity) {
    return h_n2_DyDphi->GetYaxis()->GetBinCenter(configuration->getDeltaPhiIndex(particle1, particle2) + 1);
  } else {
    return h_n2_DetaDphi->GetYaxis()->GetBinCenter(configuration->getDeltaPhiIndex(particle1, particle2) + 1);
  }
}

template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
inline int ParticlePairDiffHistos::getGlobalDeltaYEtaDeltaPhiBin(ParticleType1& p1, ParticleType2& p2)
{
  if constexpr (r == AnalysisConfiguration::kRapidity) {
    return h_n2_DyDphi->GetBin(configuration->getDeltaYEtaIndex<r>(p1, p2) + 1, configuration->getDeltaPhiIndex(p1, p2) + 1);
  } else {
    return h_n2_DetaDphi->GetBin(configuration->getDeltaYEtaIndex<r>(p1, p2) + 1, configuration->getDeltaPhiIndex(p1, p2) + 1);
  }
}

template <typename ParticleType1, typename ParticleType2>
inline float ParticlePairDiffHistos::getInvMass(ParticleType1& particle1, ParticleType2& particle2)
{
  double p1[4];
  double p2[4];
  double g[4] = {1.0, -1.0, -1.0, -1.0};
  double s = 0.0;
  float invMass = 0.0;

  particle1.getEPxPyPz(p1);
  particle2.getEPxPyPz(p2);
  for (int i = 0; i < 4; ++i) {
    s += g[i] * (p1[i] + p2[i]) * (p1[i] + p2[i]);
  }
  if (s>0){
   invMass = sqrt(s);
  }
  return invMass;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r, typename ParticleType1, typename ParticleType2>
void ParticlePairDiffHistos::getPratt(ParticleType1& particle1, ParticleType2& particle2, double& kT, double& qinv, double& qlong, double& qside, double& qout)
{
  // https://arxiv.org/abs/nucl-th/0303025
  double p1[4];
  double p2[4];
  double p[4];
  double q[4];
  double g[4] = {1.0, -1.0, -1.0, -1.0};
  double m1sq = 0.0;
  double m2sq = 0.0;
  double s = 0.0;
  double qinvsq = 0.0;

  /* we always compensate for limited acceptance */
  float rapidityCompensation = configuration->getDeltaYEtaCompensation<r>(particle1, particle2);
  float ptCompensation = configuration->getDeltaPtCompensation(particle1, particle2);

  particle1.getEPxPyPz(p1);
  particle2.getEPxPyPz(p2);

  for (int i = 0; i < 4; ++i) {
    p[i] = p1[i] + p2[i];
    q[i] = p1[i] - p2[i];
    s += g[i] * p[i] * p[i];
    m1sq += g[i] * p1[i] * p1[i];
    m2sq += g[i] * p2[i] * p2[i];
  }
  for (int i = 0; i < 4; ++i) {
    q[i] -= p[i] * (m1sq - m2sq) / s;
    qinvsq -= g[i] * q[i] * q[i];
  }

  double pT = sqrt(p[1] * p[1] + p[2] * p[2]);
  double Mlong = sqrt(s + pT * pT);

  kT = pT / 2.0;
  qinv = sqrt(qinvsq);
  if (pT > 0) {
    qlong = (p[0] * q[3] - p[3] * q[0]) / Mlong;
    qside = (p[1] * q[2] - p[2] * q[1]) / pT;
    qout = (sqrt(s) / Mlong) * (p[1] * q[1] + p[2] * q[2]) / pT;
  } else {
    qlong = (p[0] * q[3] - p[3] * q[0]) / Mlong;
    qside = sqrt(qinvsq - qlong * qlong);
    qout = 0.0;
  }
  /* we always compensate for limited acceptance */
  qinv *= rapidityCompensation * ptCompensation;
  qside *= rapidityCompensation * ptCompensation;
  qout *= rapidityCompensation * ptCompensation;
}

template <AnalysisConfiguration::RapidityPseudoRapidity r, AnalysisConfiguration::FillPairOptions options, typename ParticleType1, typename ParticleType2>
void ParticlePairDiffHistos::fill(ParticleType1& particle1, ParticleType2& particle2, double weight1, double weight2, double pTavg1, double pTavg2)
{
  int globalyetabinno = getGlobalDeltaYEtaDeltaPhiBin<r>(particle1, particle2);
  float deltayeta = getDeltaYEta<r>(particle1, particle2);
  float deltaphi = getDeltaPhi<r>(particle1, particle2);

  constexpr float pi0Mass = 0.134977;   // pi0 invariant mass
  constexpr float sidePi0MassR = 0.140; //  right mass
  constexpr float sidePi0MassL = 0.130; //  left mass
  constexpr float deltaMass = 0.002;    //  width of the mass window

  constexpr float etaMass = 0.547862;   // eta invariant mass
  constexpr float sideEtaMassR = 0.552; //  right mass
  constexpr float sideEtaMassL = 0.542; //  left mass

  constexpr float phiMass = 1.02;       // phi invariant mass
  constexpr float sidePhiMassR = 1.035; //  right mass
  constexpr float sidePhiMassL = 1.005; //  left mass

  if constexpr (r == AnalysisConfiguration::kRapidity) {
    h_n2_ptPt->Fill(particle1.pt, particle2.pt, weight1 * weight2);
    p_n2_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
    h_n2_DyDphi->AddBinContent(globalyetabinno, weight1 * weight2);
    h_ptpt_DyDphi->AddBinContent(globalyetabinno, weight1 * particle1.pt * weight2 * particle2.pt);
    h_dptdpt_DyDphi->AddBinContent(globalyetabinno, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_n2_DyDphi->SetEntries(h_n2_ptPt->GetEntries());
    h_ptpt_DyDphi->SetEntries(h_n2_ptPt->GetEntries());
    h_dptdpt_DyDphi->SetEntries(h_n2_ptPt->GetEntries());

    if constexpr ((options & AnalysisConfiguration::kFillInvariantMass) == AnalysisConfiguration::kFillInvariantMass) {
      float invMass = getInvMass(particle1, particle2);
      h_invMass->Fill(invMass);

      if (std::abs(invMass - pi0Mass) < deltaMass) {
        h_Pi0GG_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePi0MassR) < deltaMass) {
        h_Pi0GGSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePi0MassL) < deltaMass) {
        h_Pi0GGSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - etaMass) < deltaMass) {
        h_EtaGG_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sideEtaMassR) < deltaMass) {
        h_EtaGGSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sideEtaMassL) < deltaMass) {
        h_EtaGGSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - phiMass) < 2 * deltaMass) {
        h_PhiKK_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePhiMassR) < 2 * deltaMass) {
        h_PhiKKSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePhiMassL) < 2 * deltaMass) {
        h_PhiKKSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      }
    }
  } else {
    h_n2_ptPt->Fill(particle1.pt, particle2.pt, weight1 * weight2);
    p_n2_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
    h_n2_DetaDphi->AddBinContent(globalyetabinno, weight1 * weight2);
    h_ptpt_DetaDphi->AddBinContent(globalyetabinno, weight1 * particle1.pt * weight2 * particle2.pt);
    h_dptdpt_DetaDphi->AddBinContent(globalyetabinno, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_n2_DetaDphi->SetEntries(h_n2_ptPt->GetEntries());
    h_ptpt_DetaDphi->SetEntries(h_n2_ptPt->GetEntries());
    h_dptdpt_DetaDphi->SetEntries(h_n2_ptPt->GetEntries());

    if constexpr ((options & AnalysisConfiguration::kFillInvariantMass) == AnalysisConfiguration::kFillInvariantMass) {
      float invMass = getInvMass(particle1, particle2);
      h_invMass->Fill(invMass);

      if (abs(invMass - pi0Mass) < deltaMass) {
        h_Pi0GG_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePi0MassR) < deltaMass) {
        h_Pi0GGSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePi0MassL) < deltaMass) {
        h_Pi0GGSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - etaMass) < deltaMass) {
        h_EtaGG_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sideEtaMassR) < deltaMass) {
        h_EtaGGSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sideEtaMassL) < deltaMass) {
        h_EtaGGSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - phiMass) < 2 * deltaMass) {
        h_PhiKK_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePhiMassR) < 2 * deltaMass) {
        h_PhiKKSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePhiMassL) < 2 * deltaMass) {
        h_PhiKKSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      }
    }
  }

  /* TODO: this has to be templated */
  if constexpr ((options & AnalysisConfiguration::kFillPratt) == AnalysisConfiguration::kFillPratt) {
    double kT, Qinv, Qlong, Qside, Qout;
    getPratt<r>(particle1, particle2, kT, Qinv, Qlong, Qside, Qout);

    h_n2_QinvKt->Fill(kT, Qinv, weight1 * weight2);
    h_n2_QlongKt->Fill(kT, Qlong, weight1 * weight2);
    h_n2_QsideKt->Fill(kT, Qside, weight1 * weight2);
    h_n2_QoutKt->Fill(kT, Qout, weight1 * weight2);
    h_n2_QlongQsideQout->Fill(Qout, Qside, Qlong, weight1 * weight2);
    h_ptpt_QinvKt->Fill(kT, Qinv, weight1 * particle1.pt * weight2 * particle2.pt);
    h_ptpt_QlongKt->Fill(kT, Qlong, weight1 * particle1.pt * weight2 * particle2.pt);
    h_ptpt_QsideKt->Fill(kT, Qside, weight1 * particle1.pt * weight2 * particle2.pt);
    h_ptpt_QoutKt->Fill(kT, Qout, weight1 * particle1.pt * weight2 * particle2.pt);
    h_ptpt_QlongQsideQout->Fill(Qout, Qside, Qlong, weight1 * particle1.pt * weight2 * particle2.pt);
    h_dptdpt_QinvKt->Fill(kT, Qinv, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_dptdpt_QlongKt->Fill(kT, Qlong, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_dptdpt_QsideKt->Fill(kT, Qside, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_dptdpt_QoutKt->Fill(kT, Qout, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_dptdpt_QlongQsideQout->Fill(Qout, Qside, Qlong, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
  }
}

template <AnalysisConfiguration::RapidityPseudoRapidity r, AnalysisConfiguration::FillPairOptions options, typename ParticleType1, typename ParticleType2>
void ParticlePairDiffHistos::fill(ParticleType1& particle1, ParticleType2& particle2, double weight1, double weight2, TH2* pTavg1Hist, TH2* pTavg2Hist)
{
  int globalyetabinno = getGlobalDeltaYEtaDeltaPhiBin<r>(particle1, particle2);
  float deltayeta = getDeltaYEta<r>(particle1, particle2);
  float deltaphi = getDeltaPhi<r>(particle1, particle2);

  constexpr float pi0Mass = 0.134977;   // pi0 invariant mass
  constexpr float sidePi0MassR = 0.140; //  right mass
  constexpr float sidePi0MassL = 0.130; //  left mass
  constexpr float deltaMass = 0.002;    //  width of the mass window

  constexpr float etaMass = 0.547862;   // eta invariant mass
  constexpr float sideEtaMassR = 0.552; //  right mass
  constexpr float sideEtaMassL = 0.542; //  left mass

  constexpr float phiMass = 1.02;       // phi invariant mass
  constexpr float sidePhiMassR = 1.035; //  right mass
  constexpr float sidePhiMassL = 1.005; //  left mass

  float pTavg1 = pTavg1Hist->GetBinContent(pTavg1Hist->FindFixBin(particle1.eta, particle1.phi));
  float pTavg2 = pTavg2Hist->GetBinContent(pTavg2Hist->FindFixBin(particle2.eta, particle2.phi));

  if constexpr (r == AnalysisConfiguration::kRapidity) {
    h_n2_ptPt->Fill(particle1.pt, particle2.pt, weight1 * weight2);
    p_n2_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
    h_n2_DyDphi->AddBinContent(globalyetabinno, weight1 * weight2);
    h_ptpt_DyDphi->AddBinContent(globalyetabinno, weight1 * particle1.pt * weight2 * particle2.pt);
    h_dptdpt_DyDphi->AddBinContent(globalyetabinno, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_n2_DyDphi->SetEntries(h_n2_ptPt->GetEntries());
    h_ptpt_DyDphi->SetEntries(h_n2_ptPt->GetEntries());
    h_dptdpt_DyDphi->SetEntries(h_n2_ptPt->GetEntries());

    if constexpr ((options & AnalysisConfiguration::kFillInvariantMass) == AnalysisConfiguration::kFillInvariantMass) {
      float invMass = getInvMass(particle1, particle2);
      h_invMass->Fill(invMass);

      if (std::abs(invMass - pi0Mass) < deltaMass) {
        h_Pi0GG_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePi0MassR) < deltaMass) {
        h_Pi0GGSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePi0MassL) < deltaMass) {
        h_Pi0GGSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - etaMass) < deltaMass) {
        h_EtaGG_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sideEtaMassR) < deltaMass) {
        h_EtaGGSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sideEtaMassL) < deltaMass) {
        h_EtaGGSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - phiMass) < 2 * deltaMass) {
        h_PhiKK_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePhiMassR) < 2 * deltaMass) {
        h_PhiKKSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePhiMassL) < 2 * deltaMass) {
        h_PhiKKSide_DyDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      }
    }
  } else {
    h_n2_ptPt->Fill(particle1.pt, particle2.pt, weight1 * weight2);
    p_n2_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
    h_n2_DetaDphi->AddBinContent(globalyetabinno, weight1 * weight2);
    h_ptpt_DetaDphi->AddBinContent(globalyetabinno, weight1 * particle1.pt * weight2 * particle2.pt);
    h_dptdpt_DetaDphi->AddBinContent(globalyetabinno, (particle1.pt - pTavg1) * (particle2.pt - pTavg2) * weight1 * weight2);
    h_n2_DetaDphi->SetEntries(h_n2_ptPt->GetEntries());
    h_ptpt_DetaDphi->SetEntries(h_n2_ptPt->GetEntries());
    h_dptdpt_DetaDphi->SetEntries(h_n2_ptPt->GetEntries());

    if constexpr ((options & AnalysisConfiguration::kFillInvariantMass) == AnalysisConfiguration::kFillInvariantMass) {
      float invMass = getInvMass(particle1, particle2);
      h_invMass->Fill(invMass);

      if (abs(invMass - pi0Mass) < deltaMass) {
        h_Pi0GG_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePi0MassR) < deltaMass) {
        h_Pi0GGSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePi0MassL) < deltaMass) {
        h_Pi0GGSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - etaMass) < deltaMass) {
        h_EtaGG_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sideEtaMassR) < deltaMass) {
        h_EtaGGSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sideEtaMassL) < deltaMass) {
        h_EtaGGSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - phiMass) < 2 * deltaMass) {
        h_PhiKK_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePhiMassR) < 2 * deltaMass) {
        h_PhiKKSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      } else if (std::abs(invMass - sidePhiMassL) < 2 * deltaMass) {
        h_PhiKKSide_DetaDphi->Fill(deltayeta, deltaphi, weight1 * weight2);
      }
    }
  }

  /* TODO: this has to be templated */
  if constexpr ((options & AnalysisConfiguration::kFillPratt) == AnalysisConfiguration::kFillPratt) {
    double kT, Qinv, Qlong, Qside, Qout;
    getPratt<r>(particle1, particle2, kT, Qinv, Qlong, Qside, Qout);

    h_n2_QinvKt->Fill(kT, Qinv, weight1 * weight2);
    h_n2_QlongKt->Fill(kT, Qlong, weight1 * weight2);
    h_n2_QsideKt->Fill(kT, Qside, weight1 * weight2);
    h_n2_QoutKt->Fill(kT, Qout, weight1 * weight2);
    h_n2_QlongQsideQout->Fill(Qout, Qside, Qlong, weight1 * weight2);
    h_ptpt_QinvKt->Fill(kT, Qinv, weight1 * particle1.pt * weight2 * particle2.pt);
    h_ptpt_QlongKt->Fill(kT, Qlong, weight1 * particle1.pt * weight2 * particle2.pt);
    h_ptpt_QsideKt->Fill(kT, Qside, weight1 * particle1.pt * weight2 * particle2.pt);
    h_ptpt_QoutKt->Fill(kT, Qout, weight1 * particle1.pt * weight2 * particle2.pt);
    h_ptpt_QlongQsideQout->Fill(Qout, Qside, Qlong, weight1 * particle1.pt * weight2 * particle2.pt);
    h_dptdpt_QinvKt->Fill(kT, Qinv, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_dptdpt_QlongKt->Fill(kT, Qlong, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_dptdpt_QsideKt->Fill(kT, Qside, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_dptdpt_QoutKt->Fill(kT, Qout, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
    h_dptdpt_QlongQsideQout->Fill(Qout, Qside, Qlong, (weight1 * particle1.pt - pTavg1) * (weight2 * particle2.pt - pTavg2));
  }
}

#endif /* ParticlePairDiffHistos  */
