// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef WAC_ParticleHistos
#define WAC_ParticleHistos
#include "Histograms.hpp"
#include "TLorentzVector.h"
#include "AnalysisConfiguration.hpp"

class ParticleHistos : public Histograms
{
 public:
  ParticleHistos(const TString& collectionName,
                 AnalysisConfiguration* analysisConfiguration,
                 LogLevel debugLevel);
  ParticleHistos(TDirectory* dir,
                 const TString& collectionName,
                 AnalysisConfiguration* analysisConfiguration,
                 LogLevel debugLevel);
  virtual ~ParticleHistos();
  void createHistograms();
  void loadHistograms(TDirectory* dir);
  template <AnalysisConfiguration::RapidityPseudoRapidity r>
  void fill(Particle& particle, double weight);
  template <AnalysisConfiguration::RapidityPseudoRapidity r>
  void fill(TLorentzVector& p, double weight);
  void completeFill();
  void fillEventWiseInfo(float multiplicity, double nAccepted, float weight);
  void calculateAverages();

  ////////////////////////////////////////////////////////////////////////////
  // Data Members - Histograms
  ////////////////////////////////////////////////////////////////////////////
  TH1* h_n1;
  TProfile* p_n1_vsC;
  TH1* h_n1_pid;
  TH1* h_n1_pt;
  TH1* h_n1_ptXS; // 1/pt dN/dptdy
  TH1* h_n1_eta;
  TH1* h_n1_phi;
  TH2* h_n1_ptEta;
  TH2* h_n1_phiEta;
  TH3* h_n1_ptPhiEta;
  TH1* h_spt_phi;
  TH1* h_spt_eta;
  TH2* h_spt_phiEta;
  TH1* h_pt_phi;
  TH1* h_pt_eta;
  TH2* h_pt_phiEta;

  TH1* h_n1_y;
  TH2* h_n1_ptY;
  TH2* h_n1_phiY;
  TH3* h_n1_ptPhiY;
  TH1* h_spt_y;
  TH2* h_spt_phiY;
  TH1* h_pt_y;
  TH2* h_pt_phiY;

  ClassDef(ParticleHistos, 4)
};

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void ParticleHistos::fill(Particle& particle, double weight)
{
  double pt = particle.pt;
  double phi = particle.phi;
  if (phi < 0)
    phi += TMath::TwoPi();

  h_n1_pid->Fill(TString::Format("%ld", particle.pid), 1);
  h_n1_pt->Fill(pt, weight);
  h_n1_ptXS->Fill(pt, weight / pt);
  if constexpr (r == AnalysisConfiguration::kRapidity) {
    double y = particle.y;
    // delayed h_n1_y->Fill(y, weight);
    h_n1_ptY->Fill(y, pt, weight);
    h_n1_phiY->Fill(y, phi, weight);
    // delatyed h_spt_y->Fill(y, pt * weight);
    h_spt_phiY->Fill(y, phi, pt * weight);
  } else {
    double eta = particle.eta;
    // delayed fill h_n1_eta    ->Fill(eta, weight);
    // delayed fill h_n1_phi    ->Fill(phi, weight);
    h_n1_ptEta->Fill(eta, pt, weight);
    h_n1_phiEta->Fill(eta, phi, weight);
    // delayed fill h_spt_phi    ->Fill(phi, pt*weight);
    // delayed fill h_spt_eta    ->Fill(eta, pt*weight);
    h_spt_phiEta->Fill(eta, phi, pt * weight);
  }
}

template <AnalysisConfiguration::RapidityPseudoRapidity r>
void ParticleHistos::fill(TLorentzVector& p, double weight)
{
  double pt = p.Pt();
  double phi = p.Phi();
  if (phi < 0)
    phi += TMath::TwoPi();

  h_n1_pt->Fill(pt, weight);
  h_n1_ptXS->Fill(pt, weight / pt);

  if constexpr (r == AnalysisConfiguration::kRapidity) {
    double y = p.Rapidity();
    // delayed fill h_n1_y      ->Fill(y, weight);
    h_n1_ptY->Fill(y, pt, weight);
    h_n1_phiY->Fill(y, phi, weight);
    // delayed fill h_spt_y     ->Fill(y, pt*weight);
    h_spt_phiY->Fill(y, phi, pt * weight);
  } else {
    double eta = p.Eta();
    // delayed fill h_n1_eta     ->Fill(eta, weight);
    // delayed fill h_n1_phi     ->Fill(phi, weight);
    h_n1_ptEta->Fill(eta, pt, weight);
    h_n1_phiEta->Fill(eta, phi, weight);
    // delayed fill h_spt_phi    ->Fill(phi, pt*weight);
    // delayed fill h_spt_eta    ->Fill(eta, pt*weight);
    h_spt_phiEta->Fill(eta, phi, pt * weight);
  }
}

#endif /* WAC_ParticleHistos  */
