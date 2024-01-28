//
//  ParticleHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include "TLorentzVector.h"
#include "ParticleHistos.hpp"
ClassImp(ParticleHistos)

  ParticleHistos::ParticleHistos(const TString& name,
                                 AnalysisConfiguration* configuration,
                                 LogLevel debugLevel)
  : Histograms(name, configuration, 100, debugLevel)
{
  initialize();
}

ParticleHistos::ParticleHistos(TDirectory* dir,
                               const TString& name,
                               AnalysisConfiguration* configuration,
                               LogLevel debugLevel)
  : Histograms(name, configuration, 100, debugLevel)
{
  loadHistograms(dir);
}

ParticleHistos::~ParticleHistos()
{
  // deleteHistograms();
}

// for now use the same boundaries for eta and y histogram
void ParticleHistos::createHistograms()
{
  AnalysisConfiguration& ac = *getConfiguration();
  TString bn = getHistoBaseName();
  h_n1 = createHistogram(bn + TString("n1"), 1000, -0.5, 999.5, "n_1", "N", scaled, saved, plotted, notPrinted);
  p_n1_vsC = createProfile(bn + TString("n1_vsC"), 101, -0.5,100.5, "Centrality/Multiplicity (%)", "#LTn_{1}#GT", saved, notPlotted, notPrinted);
  h_n1_pid = createHistogram(bn + TString("n1_pid"), 100, 0.5, 100.5, "PID", "N", scaled, saved, plotted, notPrinted);
  if (ac.bin_edges_pt.size() > 0) {
    h_n1_pt = createHistogram(bn + TString("n1_pt"), ac.bin_edges_pt, "p_{T}", "N", scaled, saved, plotted, notPrinted);
    h_n1_ptXS = createHistogram(bn + TString("n1_ptXS"), ac.bin_edges_pt, "p_{T}", "1/p_{T} dN/p_{T}", scaled, saved, plotted, notPrinted);
  } else {
    h_n1_pt = createHistogram(bn + TString("n1_pt"), ac.nBins_pt, ac.min_pt, ac.max_pt, "p_{T}", "N", scaled, saved, plotted, notPrinted);
    h_n1_ptXS = createHistogram(bn + TString("n1_ptXS"), ac.nBins_pt, ac.min_pt, ac.max_pt, "p_{T}", "1/p_{T} dN/p_{T}", scaled, saved, plotted, notPrinted);
  }
  h_n1_phi = createHistogram(bn + TString("n1_phi"), ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi", "N", scaled, saved, plotted, notPrinted);
  h_spt_phi = createHistogram(bn + TString("sumpt1_phi"), ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi", "#sum p_{T}", scaled, saved, plotted, notPrinted);
  h_pt_phi = createHistogram(bn + TString("pt1_phi"), ac.nBins_phi, ac.min_phi, ac.max_phi, "#varphi", "#LTp_{T}#GT", notScaled, saved, plotted, notPrinted);

  if (ac.fillYorEta == ac.kPseudorapidity) {
    h_n1_eta = createHistogram(bn + TString("n1_eta"), ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta", "N", scaled, saved, plotted, notPrinted);
    if (ac.bin_edges_pt.size() > 0) {
      h_n1_ptEta = createHistogram(bn + TString("n1_ptEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.bin_edges_pt, "#eta", "p_{T}", "N", scaled, saved, plotted, notPrinted);
    } else {
      h_n1_ptEta = createHistogram(bn + TString("n1_ptEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.nBins_pt, ac.min_pt, ac.max_pt, "#eta", "p_{T}", "N", scaled, saved, plotted, notPrinted);
    }
    h_n1_phiEta = createHistogram(bn + TString("n1_phiEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.nBins_phi, ac.min_phi, ac.max_phi, "#eta", "#varphi", "N", scaled, saved, plotted, notPrinted);
    h_spt_eta = createHistogram(bn + TString("sumpt1_eta"), ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta", "#sum p_{T}", scaled, saved, plotted, notPrinted);
    h_spt_phiEta = createHistogram(bn + TString("sumpt1_phiEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.nBins_phi, ac.min_phi, ac.max_phi, "#eta", "#varphi", "#sum p_{T}", scaled, saved, notPlotted, notPrinted);
    h_pt_eta = createHistogram(bn + TString("pt1_eta"), ac.nBins_eta, ac.min_eta, ac.max_eta, "#eta", "#LTp_{T}#GT", notScaled, saved, plotted, notPrinted);
    h_pt_phiEta = createHistogram(bn + TString("pt1_phiEta"), ac.nBins_eta, ac.min_eta, ac.max_eta, ac.nBins_phi, ac.min_phi, ac.max_phi, "#eta", "#varphi", "#LTp_{T}#GT", notScaled, saved, notPlotted, notPrinted);
  }
  if (ac.fillYorEta == ac.kRapidity) {
    h_n1_y = createHistogram(bn + TString("n1_y"), ac.nBins_y, ac.min_y, ac.max_y, "y", "N", scaled, saved, plotted, notPrinted);
    if (ac.bin_edges_pt.size() > 0) {
      h_n1_ptY = createHistogram(bn + TString("n1_ptY"), ac.nBins_y, ac.min_y, ac.max_y, ac.bin_edges_pt, "#eta", "p_{T}", "N", scaled, saved, notPlotted, notPrinted);
    } else {
      h_n1_ptY = createHistogram(bn + TString("n1_ptY"), ac.nBins_y, ac.min_y, ac.max_y, ac.nBins_pt, ac.min_pt, ac.max_pt, "#eta", "p_{T}", "N", scaled, saved, notPlotted, notPrinted);
    }
    h_n1_phiY = createHistogram(bn + TString("n1_phiY"), ac.nBins_y, ac.min_y, ac.max_y, ac.nBins_phi, ac.min_phi, ac.max_phi, "y", "#varphi", "N", scaled, saved, notPlotted, notPrinted);
    h_spt_y = createHistogram(bn + TString("sumpt1_y"), ac.nBins_y, ac.min_y, ac.max_y, "y", "#sum p_{T}", scaled, saved, notPlotted, notPrinted);
    h_spt_phiY = createHistogram(bn + TString("sumpt1_phiY"), ac.nBins_y, ac.min_y, ac.max_y, ac.nBins_phi, ac.min_phi, ac.max_phi, "y", "#varphi", "#sum p_{T}", scaled, saved, notPlotted, notPrinted);
    h_pt_y = createHistogram(bn + TString("pt1_y"), ac.nBins_y, ac.min_y, ac.max_y, "y", "#LTp_{T}#GT", notScaled, saved, notPlotted, notPrinted);
    h_pt_phiY = createHistogram(bn + TString("pt1_phiY"), ac.nBins_y, ac.min_y, ac.max_y, ac.nBins_phi, ac.min_phi, ac.max_phi, "y", "#varphi", "#LTp_{T}#GT", notScaled, saved, notPlotted, notPrinted);
  }

  if (ac.fill3D) {
    if (ac.fillYorEta == ac.kPseudorapidity) {
      h_n1_ptPhiEta = createHistogram(bn + TString("n1_ptPhiEta"),
                                      ac.nBins_eta, ac.min_eta, ac.max_eta,
                                      ac.nBins_phi, ac.min_phi, ac.max_phi,
                                      ac.nBins_pt, ac.min_pt, ac.max_pt,
                                      "#eta", "#varphi", "p_{T}", "N", scaled, saved, notPlotted, notPrinted);
    }
    if (ac.fillYorEta == ac.kRapidity) {
      h_n1_ptPhiY = createHistogram(bn + TString("n1_ptPhiY"),
                                    100, ac.min_y, ac.max_y,
                                    100, ac.min_phi, ac.max_phi,
                                    100, ac.min_pt, ac.max_pt,
                                    "y", "#varphi", "p_{T}", "N", scaled, saved, notPlotted, notPrinted);
    }
  }
}

//________________________________________________________________________
void ParticleHistos::loadHistograms(TDirectory* dir)
{
  if (!dir) {
    if (reportFatal())
      cout << "-Fatal- Attempting to load ParticleHistos from an invalid directory pointer" << endl;
    return;
  }
  dir->cd();

  AnalysisConfiguration& ac = *getConfiguration();
  TString bn = getHistoBaseName();
  h_n1 = loadH1(dir, bn + TString("n1"), true);
  p_n1_vsC = loadProfile(dir, bn + TString("n1_vsC"), false);
  h_n1_pid = loadH1(dir, bn + TString("n1_pid"), true);
  h_n1_pt = loadH1(dir, bn + TString("n1_pt"), true);
  h_n1_ptXS = loadH1(dir, bn + TString("n1_ptXS"), true);
  h_n1_phi = loadH1(dir, bn + TString("n1_phi"), true);
  h_spt_phi = loadH1(dir, bn + TString("sumpt1_phi"), true);
  h_pt_phi = loadH1(dir, bn + TString("pt1_phi"), false, true);
  if (ac.fillYorEta == ac.kPseudorapidity) {
    h_n1_eta = loadH1(dir, bn + TString("n1_eta"), true);
    h_n1_ptEta = loadH2(dir, bn + TString("n1_ptEta"), true);
    h_n1_phiEta = loadH2(dir, bn + TString("n1_phiEta"), true);
    h_spt_phiEta = loadH2(dir, bn + TString("sumpt1_phiEta"), true);
    h_spt_eta = loadH1(dir, bn + TString("sumpt1_eta"), true);
    h_pt_phiEta = loadH2(dir, bn + TString("pt1_phiEta"), false, true);
    h_pt_eta = loadH1(dir, bn + TString("pt1_eta"), false, true);
    if (ac.fill3D) {
      h_n1_ptPhiEta = loadH3(dir, bn + TString("n1_ptPhiEta"), true);
    }
  }
  if (ac.fillYorEta == ac.kRapidity) {
    h_n1_y = loadH1(dir, bn + TString("n1_y"), true);
    h_n1_ptY = loadH2(dir, bn + TString("n1_ptY"), true);
    h_n1_phiY = loadH2(dir, bn + TString("n1_phiY"), true);
    h_spt_phiY = loadH2(dir, bn + TString("sumpt1_phiY"), true);
    h_spt_y = loadH1(dir, bn + TString("sumpt1_y"), true);
    h_pt_phiY = loadH2(dir, bn + TString("pt1_phiY"), false, true);
    h_pt_y = loadH1(dir, bn + TString("pt1_y"), false, true);
    if (ac.fill3D) {
      h_n1_ptPhiY = loadH3(dir, bn + TString("n1_ptPhiY"), true);
    }
  }
  /* the histograms are not owned */
  bOwnTheHistograms = false;
  return;
}

void ParticleHistos::fillEventWiseInfo(float multiplicity, double nAccepted, float weight)
{
  h_n1->Fill(nAccepted, weight);
  p_n1_vsC->Fill(multiplicity, nAccepted, weight);
}

// complete filling the addicional histograms by projecting the
// higher dimensional ones
void ParticleHistos::completeFill()
{
  AnalysisConfiguration& ac = *getConfiguration();

  if (ac.fillYorEta == ac.kPseudorapidity) {
    TH1* h_eta = h_n1_phiEta->ProjectionX();
    h_n1_eta->Reset();
    h_n1_eta->Add(h_eta);
    TH1* h_phi = h_n1_phiEta->ProjectionY();
    h_n1_phi->Reset();
    h_n1_phi->Add(h_phi);
    TH1* h_pteta = h_spt_phiEta->ProjectionX();
    h_spt_eta->Reset();
    h_spt_eta->Add(h_pteta);
    TH1* h_ptphi = h_spt_phiEta->ProjectionY();
    h_spt_phi->Reset();
    h_spt_phi->Add(h_ptphi);
    delete h_ptphi;
    delete h_pteta;
    delete h_phi;
    delete h_eta;
  }

  if (ac.fillYorEta == ac.kRapidity) {
    TH1* h_y = h_n1_phiY->ProjectionX();
    h_n1_y->Reset();
    h_n1_y->Add(h_y);
    TH1* h_phi = h_n1_phiY->ProjectionY();
    h_n1_phi->Reset();
    h_n1_phi->Add(h_phi);
    TH1* h_pty = h_spt_phiY->ProjectionX();
    h_spt_y->Reset();
    h_spt_y->Add(h_pty);
    TH1* h_ptphi = h_spt_phiY->ProjectionY();
    h_spt_phi->Reset();
    h_spt_phi->Add(h_ptphi);
    delete h_ptphi;
    delete h_pty;
    delete h_phi;
    delete h_y;
  }
}

void ParticleHistos::calculateAverages()
{
  AnalysisConfiguration& ac = *getConfiguration();
  calculateAveragePt(h_spt_phi, h_n1_phi, h_pt_phi);
  if (ac.fillYorEta == ac.kRapidity) {
    calculateAveragePt(h_spt_y, h_n1_y, h_pt_y);
    calculateAveragePtH2(h_spt_phiY, h_n1_phiY, h_pt_phiY);
  } else {
    calculateAveragePt(h_spt_eta, h_n1_eta, h_pt_eta);
    calculateAveragePtH2(h_spt_phiEta, h_n1_phiEta, h_pt_phiEta);
  }
}
