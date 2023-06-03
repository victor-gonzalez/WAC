//
//  DerivedParticleHistos.cpp
//  MyMC
//
//  Created by Claude Pruneau on 9/23/16.
//  Copyright © 2016 Claude Pruneau. All rights reserved.
//

#include "ParticlePairDerivedDiffHistos.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////////////////////////////////////////////
ParticlePairDerivedDiffHistos::ParticlePairDerivedDiffHistos(const TString& name,
                                                             AnalysisConfiguration* configuration,
                                                             LogLevel debugLevel)
  : ParticlePairDiffHistos(name, configuration, debugLevel),
    h_n1n1_ptPt(nullptr),
    h_n1n1_DetaDphi(nullptr),
    h_pt1pt1_DetaDphi(nullptr),
    h_R2_DetaDphi(nullptr),
    h_P2_DetaDphi(nullptr),
    h_G2_DetaDphi(nullptr),
    h_R2_DetaDphi_shft(nullptr),
    h_P2_DetaDphi_shft(nullptr),
    h_G2_DetaDphi_shft(nullptr),
    h_R2bf12_DetaDphi_shft(nullptr),
    h_PrattBf_DetaDphi_shft(nullptr),
    p_PrattBf_DetaDphi_shft(nullptr),
    h_n1n1_DyDphi(nullptr),
    h_pt1pt1_DyDphi(nullptr),
    h_R2_DyDphi(nullptr),
    h_P2_DyDphi(nullptr),
    h_G2_DyDphi(nullptr),
    h_R2_DyDphi_shft(nullptr),
    h_P2_DyDphi_shft(nullptr),
    h_G2_DyDphi_shft(nullptr),
    h_R2bf12_DyDphi_shft(nullptr),
    h_PrattBf_DyDphi_shft(nullptr),
    p_PrattBf_DyDphi_shft(nullptr)
{
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::CTOR() Started." << endl;
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::CTOR() Completed." << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// CTOR
//////////////////////////////////////////////////////////////////////////////////////////////////////
ParticlePairDerivedDiffHistos::ParticlePairDerivedDiffHistos(TDirectory* dir,
                                                             const TString& name,
                                                             AnalysisConfiguration* configuration,
                                                             LogLevel debugLevel)
  : ParticlePairDiffHistos(dir, name, configuration, debugLevel),
    h_n1n1_ptPt(nullptr),
    h_n1n1_DetaDphi(nullptr),
    h_pt1pt1_DetaDphi(nullptr),
    h_R2_DetaDphi(nullptr),
    h_P2_DetaDphi(nullptr),
    h_G2_DetaDphi(nullptr),
    h_R2_DetaDphi_shft(nullptr),
    h_P2_DetaDphi_shft(nullptr),
    h_G2_DetaDphi_shft(nullptr),
    h_R2bf12_DetaDphi_shft(nullptr),
    h_PrattBf_DetaDphi_shft(nullptr),
    p_PrattBf_DetaDphi_shft(nullptr),
    h_n1n1_DyDphi(nullptr),
    h_pt1pt1_DyDphi(nullptr),
    h_R2_DyDphi(nullptr),
    h_P2_DyDphi(nullptr),
    h_G2_DyDphi(nullptr),
    h_R2_DyDphi_shft(nullptr),
    h_P2_DyDphi_shft(nullptr),
    h_G2_DyDphi_shft(nullptr),
    h_R2bf12_DyDphi_shft(nullptr),
    h_PrattBf_DyDphi_shft(nullptr),
    p_PrattBf_DyDphi_shft(nullptr)
{
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::CTOR(TFile * inputFile,...) Started." << endl;
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::CTOR(TFile * inputFile,...) Completed." << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// DTOR
//////////////////////////////////////////////////////////////////////////////////////////////////////
ParticlePairDerivedDiffHistos::~ParticlePairDerivedDiffHistos()
{
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::DTOR() No ops." << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialize the histograms
//////////////////////////////////////////////////////////////////////////////////////////////////////
void ParticlePairDerivedDiffHistos::createDerivedHistograms()
{
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::initialize() starting" << endl;
  TString name;
  AnalysisConfiguration& ac = *(AnalysisConfiguration*)getConfiguration();
  TString bn = getHistoBaseName();

  if (ac.bin_edges_pt.size() > 0) {
    h_n1n1_ptPt = createHistogram(bn + TString("n1n1_ptPt"), ac.bin_edges_pt, ac.bin_edges_pt, "p_{T,1}", "p_{T,2}", "<n_{1}><n_{2}>", notScaled, saved, plotted, notPrinted);
  } else {
    h_n1n1_ptPt = createHistogram(bn + TString("n1n1_ptPt"), ac.nBins_pt, ac.min_pt, ac.max_pt, ac.nBins_pt, ac.min_pt, ac.max_pt, "p_{T,1}", "p_{T,2}", "<n_{1}><n_{2}>", notScaled, saved, plotted, notPrinted);
  }

  if (ac.fillYorEta == ac.kPseudorapidity) {
    h_n1n1_DetaDphi = createHistogram(bn + TString("n1n1_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "<n_{1}><n_{1}>", notScaled, saved, plotted, notPrinted);
    h_pt1pt1_DetaDphi = createHistogram(bn + TString("pt1pt1_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "pt1pt1", notScaled, saved, plotted, notPrinted);

    h_R2_DetaDphi = createHistogram(bn + TString("R2_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "R_{2}", notScaled, saved, plotted, notPrinted);
    h_P2_DetaDphi = createHistogram(bn + TString("P2_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "P_{2}", notScaled, saved, plotted, notPrinted);
    h_G2_DetaDphi = createHistogram(bn + TString("G2_DetaDphi"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "G_{2}", notScaled, saved, plotted, notPrinted);

    h_R2_DetaDphi_shft = createHistogram(bn + TString("R2_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "R_{2}", notScaled, saved, plotted, notPrinted);
    h_P2_DetaDphi_shft = createHistogram(bn + TString("P2_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "P_{2}", notScaled, saved, plotted, notPrinted);
    h_G2_DetaDphi_shft = createHistogram(bn + TString("G2_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "G_{2}", notScaled, saved, plotted, notPrinted);

    h_R2bf12_DetaDphi_shft = createHistogram(bn + TString("R2bf12_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "#rho_{1_{1}}R_{2}^{12}", notScaled, saved, plotted, notPrinted);
    h_PrattBf_DetaDphi_shft = createHistogram(bn + TString("N2PratBf_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta#eta", "#Delta#varphi", "#rho_{2}^{12}/#rho_{1_{2}}", notScaled, saved, plotted, notPrinted);
    p_PrattBf_DetaDphi_shft = createProfile(bn + TString("PrattBf_DetaDphi_shft"), ac.nBins_Deta, ac.min_Deta, ac.max_Deta, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta#eta", "#Delta#varphi", "#rho_{2}^{12}/#rho_{1_{2}}", notScaled, saved, plotted, notPrinted);
  }

  if (ac.fillYorEta == ac.kRapidity) {
    h_n1n1_DyDphi = createHistogram(bn + TString("n1n1_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "<n_{1}><n_{1}>", notScaled, saved, plotted, notPrinted);
    h_pt1pt1_DyDphi = createHistogram(bn + TString("pt1pt1_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "pt1pt1", notScaled, saved, notPlotted, notPrinted);

    h_R2_DyDphi = createHistogram(bn + TString("R2_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "R_{2}", notScaled, saved, notPlotted, notPrinted);
    h_P2_DyDphi = createHistogram(bn + TString("P2_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "P_{2}", notScaled, saved, notPlotted, notPrinted);
    h_G2_DyDphi = createHistogram(bn + TString("G2_DyDphi"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "G_{2}", notScaled, saved, notPlotted, notPrinted);

    h_R2_DyDphi_shft = createHistogram(bn + TString("R2_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "R_{2}", notScaled, saved, plotted, notPrinted);
    h_P2_DyDphi_shft = createHistogram(bn + TString("P2_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "P_{2}", notScaled, saved, notPlotted, notPrinted);
    h_G2_DyDphi_shft = createHistogram(bn + TString("G2_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "G_{2}", notScaled, saved, notPlotted, notPrinted);

    h_R2bf12_DyDphi_shft = createHistogram(bn + TString("R2bf12_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "#rho_{1_{1}}R_{2}^{12}", notScaled, saved, notPlotted, notPrinted);
    h_PrattBf_DyDphi_shft = createHistogram(bn + TString("N2PrattBf_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi_shft, ac.max_Dphi_shft, "#Delta y", "#Delta#varphi", "#rho_{2}^{12}/#rho_{1_{2}}", notScaled, saved, notPlotted, notPrinted);
    p_PrattBf_DyDphi_shft = createProfile(bn + TString("PrattBf_DyDphi_shft"), ac.nBins_Dy, ac.min_Dy, ac.max_Dy, ac.nBins_Dphi, ac.min_Dphi, ac.max_Dphi, "#Delta y", "#Delta#varphi", "#rho_{2}^{12}/#rho_{1_{2}}", notScaled, saved, notPlotted, notPrinted);
  }
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::initialize() completed." << endl;
}

void ParticlePairDerivedDiffHistos::reduce(const TH2* source, TH2* target, int nEtaBins, int nPhiBins)
{
  double v1, ev1;
  int dPhi, dEta, iPhi, iEta, jPhi, jEta, i, j;
  int nWrk = nPhiBins * (2 * nEtaBins - 1);
  int index;
  double* numerator = new double[nWrk];
  double* numeratorErr = new double[nWrk];
  for (int k = 0; k < nWrk; ++k) {
    numerator[k] = 0;
    numeratorErr[k] = 0;
  }

  i = 1;
  for (iEta = 0; iEta < nEtaBins; ++iEta) {
    for (iPhi = 0; iPhi < nPhiBins; ++iPhi) {
      j = 1;
      for (jEta = 0; jEta < nEtaBins; ++jEta) {
        for (jPhi = 0; jPhi < nPhiBins; ++jPhi) {
          dPhi = iPhi - jPhi;
          if (dPhi < 0)
            dPhi += nPhiBins;
          dPhi += 1;
          dEta = iEta - jEta + nEtaBins;
          v1 = source->GetBinContent(i, j);
          ev1 = source->GetBinError(i, j);
          index = (dEta - 1) * nPhiBins + dPhi - 1;
          numerator[index] += v1;
          numeratorErr[index] += ev1 * ev1;
          ++j;
        }
      }
      ++i;
    }
  }

  for (dEta = 0; dEta < 2 * nEtaBins - 1; ++dEta) {
    for (dPhi = 0; dPhi < nPhiBins; ++dPhi) {
      index = dEta * nPhiBins + dPhi;
      v1 = numerator[index];
      ev1 = numeratorErr[index];
      target->SetBinContent(dEta + 1, dPhi + 1, v1);
      target->SetBinError(dEta + 1, dPhi + 1, sqrt(ev1));
    }
  }
  delete[] numerator;
  delete[] numeratorErr;
}

void ParticlePairDerivedDiffHistos::downscale(const TH2* source, TH2* target, int nEtaBins, int nPhiBins)
{
  /* downscale to compensate the limited acceptance */
  int dPhi, dEta, iPhi, iEta, jPhi, jEta, i, j;
  int nWrk = nPhiBins * (2 * nEtaBins - 1);
  int index;
  double* denominator = new double[nWrk];
  for (int k = 0; k < nWrk; ++k) {
    denominator[k] = 0;
  }

  i = 1;
  for (iEta = 0; iEta < nEtaBins; ++iEta) {
    for (iPhi = 0; iPhi < nPhiBins; ++iPhi) {
      j = 1;
      for (jEta = 0; jEta < nEtaBins; ++jEta) {
        for (jPhi = 0; jPhi < nPhiBins; ++jPhi) {
          dPhi = iPhi - jPhi;
          if (dPhi < 0)
            dPhi += nPhiBins;
          dPhi += 1;
          dEta = iEta - jEta + nEtaBins;
          index = (dEta - 1) * nPhiBins + dPhi - 1;
          denominator[index] += 1.;
          ++j;
        }
      }
      ++i;
    }
  }

  double max = -1;
  for (int k = 0; k < nWrk; ++k) {
    if (max < denominator[k]) {
      max = denominator[k];
    }
  }
  for (int k = 0; k < nWrk; ++k) {
    denominator[k] /= max;
  }

  target->Reset();
  for (dEta = 0; dEta < 2 * nEtaBins - 1; ++dEta) {
    for (dPhi = 0; dPhi < nPhiBins; ++dPhi) {
      index = dEta * nPhiBins + dPhi;
      double v = denominator[index];
      target->SetBinContent(dEta + 1, dPhi + 1, source->GetBinContent(dEta + 1, dPhi + 1) / v);
      target->SetBinError(dEta + 1, dPhi + 1, source->GetBinError(dEta + 1, dPhi + 1) / v);
    }
  }
  delete[] denominator;
}

double ParticlePairDerivedDiffHistos::downscaleMax(int nEtaBins, int nPhiBins)
{
  /* maximum downscale to compensate the limited acceptance */
  /* always presesrve the integral                  */
  int dPhi, dEta, iPhi, iEta, jPhi, jEta, i, j;
  int nWrk = nPhiBins * (2 * nEtaBins - 1);
  int index;
  double* denominator = new double[nWrk];
  for (int k = 0; k < nWrk; ++k) {
    denominator[k] = 0;
  }

  i = 1;
  for (iEta = 0; iEta < nEtaBins; ++iEta) {
    for (iPhi = 0; iPhi < nPhiBins; ++iPhi) {
      j = 1;
      for (jEta = 0; jEta < nEtaBins; ++jEta) {
        for (jPhi = 0; jPhi < nPhiBins; ++jPhi) {
          dPhi = iPhi - jPhi;
          if (dPhi < 0)
            dPhi += nPhiBins;
          dPhi += 1;
          dEta = iEta - jEta + nEtaBins;
          index = (dEta - 1) * nPhiBins + dPhi - 1;
          denominator[index] += 1.;
          ++j;
        }
      }
      ++i;
    }
  }

  double max = -1;
  for (int k = 0; k < nWrk; ++k) {
    if (max < denominator[k]) {
      max = denominator[k];
    }
  }
  delete[] denominator;
  return max;
}

// Histograms from ParticlePairDerivedDiffHistos must be normalized "per event" before calling this function
void ParticlePairDerivedDiffHistos::calculateDerivedHistograms(ParticleHistos* part1Histos, ParticleHistos* part2Histos, double bincorrection)
{
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::calculateDerivedHistograms() started." << endl;

  if (configuration->fillYorEta == configuration->kPseudorapidity) {
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // R2 related histograms
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    TH2F* h_n1n1_phiEtaPhiEta = new TH2F((getHistoBaseName() + TString("n1n1_phiEtaPhiEta")).Data(),
                                         "<n_{1}><n_{1}>;#eta_{1} x #varphi_{1};#eta_{2} x #varphi_{2}",
                                         configuration->nBins_phiEta, 0.0, configuration->nBins_phiEta,
                                         configuration->nBins_phiEta, 0.0, configuration->nBins_phiEta);
    TH2F* h_pt1pt1_phiEtaPhiEta = new TH2F((getHistoBaseName() + TString("pt1pt1_phiEtaPhiEta")).Data(),
                                           "<pT_{1}><pT_{1}>;#eta_{1} x #varphi_{1};#eta_{2} x #varphi_{2}",
                                           configuration->nBins_phiEta, 0.0, configuration->nBins_phiEta,
                                           configuration->nBins_phiEta, 0.0, configuration->nBins_phiEta);

    calculateN1N1_H2H2H2(part1Histos->h_n1_phiEta, part2Histos->h_n1_phiEta, h_n1n1_phiEtaPhiEta, 1.0, 1.0);
    calculateN1N1_H1H1H2(part1Histos->h_n1_pt, part2Histos->h_n1_pt, h_n1n1_ptPt, 1.0, 1.0);

    calculateN1N1_H2H2H2(part1Histos->h_spt_phiEta, part2Histos->h_spt_phiEta, h_pt1pt1_phiEtaPhiEta, 1.0, 1.0);

    reduce(h_n1n1_phiEtaPhiEta, h_n1n1_DetaDphi, configuration->nBins_eta, configuration->nBins_phi);
    reduce(h_pt1pt1_phiEtaPhiEta, h_pt1pt1_DetaDphi, configuration->nBins_eta, configuration->nBins_phi);

    TF2* minusone = new TF2("minusone", "-1",
                            h_R2_DetaDphi->GetXaxis()->GetBinLowEdge(1),
                            h_R2_DetaDphi->GetXaxis()->GetBinUpEdge(h_R2_DetaDphi->GetNbinsX()),
                            h_R2_DetaDphi->GetYaxis()->GetBinLowEdge(1),
                            h_R2_DetaDphi->GetYaxis()->GetBinUpEdge(h_R2_DetaDphi->GetNbinsY()));
    double pt_1_average = part1Histos->h_pt_phiEta->Integral() / (configuration->nBins_eta * configuration->nBins_phi);
    double pt_2_average = part2Histos->h_pt_phiEta->Integral() / (configuration->nBins_eta * configuration->nBins_phi);
    if (reportInfo()) {
      cout << "calculateDerivedDiffHistograms(...) avgPt1 eta = " << pt_1_average << endl;
      cout << "calculateDerivedDiffHistograms(...) avgPt2 eta = " << pt_2_average << endl;
    }

    /* calculate R2 */
    h_R2_DetaDphi->Reset();
    h_R2_DetaDphi->Add(h_n2_DetaDphi, bincorrection);
    h_R2_DetaDphi->Divide(h_n1n1_DetaDphi);
    h_R2_DetaDphi->Add(minusone);
    shiftY(*h_R2_DetaDphi, *h_R2_DetaDphi_shft, configuration->nBins_Dphi_shft);

    /* calculate G2 */
    h_G2_DetaDphi->Reset();
    h_G2_DetaDphi->Add(h_ptpt_DetaDphi, bincorrection);
    h_G2_DetaDphi->Add(h_pt1pt1_DetaDphi, -1.0);
    h_G2_DetaDphi->Divide(h_n1n1_DetaDphi);
    h_G2_DetaDphi->Scale(1.0 / (pt_1_average * pt_2_average));
    shiftY(*h_G2_DetaDphi, *h_G2_DetaDphi_shft, configuration->nBins_Dphi_shft);

    /* calculate P2 */
    h_P2_DetaDphi->Reset();
    h_P2_DetaDphi->Add(h_dptdpt_DetaDphi, 1.0 / (pt_1_average * pt_2_average));
    h_P2_DetaDphi->Divide(h_n2_DetaDphi);
    shiftY(*h_P2_DetaDphi, *h_P2_DetaDphi_shft, configuration->nBins_Dphi_shft);

    /* calculate BF component from R2 */
    double n1_2 = part2Histos->h_n1_phiEta->Integral();
    double volume = kTWOPI * (part2Histos->h_n1_phiEta->GetXaxis()->GetBinUpEdge(part2Histos->h_n1_phiEta->GetNbinsX()) - part2Histos->h_n1_phiEta->GetXaxis()->GetBinLowEdge(1));
    h_R2bf12_DetaDphi_shft->Reset();
    h_R2bf12_DetaDphi_shft->Add(h_R2_DetaDphi_shft, n1_2 / volume);

    /* calculate Pratt's BF component from profiles */
    /* in principle the second single division should be automatically incorporated by the profiles when the same sign combination is substracted */
    /* we havo to compensate for the extra downscale by the profile mechanism and convert it to a density distribution */
    /* No, we don't! all the compensation should be taken care by the TProfile mechanism */
    p_PrattBf_DetaDphi_shft->Reset();
    p_PrattBf_DetaDphi_shft->Add(p_n2_DetaDphi, 1.0);

    /* calculate Pratt's BF component from n2 and n1 histograms */
    TH2* tbf = (TH2*)h_n2_DetaDphi->Clone(TString::Format("%s_clone", h_PrattBf_DetaDphi_shft->GetName()).Data());
    tbf->Reset();
    tbf->Add(h_n2_DetaDphi, h_n1n1_DetaDphi, 1.0, -1.0);
    shiftY(*tbf, *h_PrattBf_DetaDphi_shft, configuration->nBins_Dphi_shft);
    h_PrattBf_DetaDphi_shft->Scale(1.0 / (n1_2 * h_PrattBf_DetaDphi_shft->GetXaxis()->GetBinWidth(1) * h_PrattBf_DetaDphi_shft->GetYaxis()->GetBinWidth(1)));
    delete tbf;

    delete h_n1n1_phiEtaPhiEta;
    delete h_pt1pt1_phiEtaPhiEta;
    delete minusone;
  }

  if (configuration->fillYorEta == configuration->kRapidity) {
    TH2F* h_n1n1_phiYPhiY = new TH2F((getHistoBaseName() + TString("n1n1_phiYPhiY")).Data(),
                                     "<n_{1}><n_{1}>;y_{1} x #varphi_{1};y_{2} x #varphi_{2}",
                                     configuration->nBins_phiY, 0.0, configuration->nBins_phiY,
                                     configuration->nBins_phiY, 0.0, configuration->nBins_phiY);
    TH2F* h_pt1pt1_phiYPhiY = new TH2F((getHistoBaseName() + TString("pt1pt1_phiYPhiY")).Data(),
                                       "<pT_{1}><pT_{1}>;y_{1} x #varphi_{1};y_{2} x #varphi_{2}",
                                       configuration->nBins_phiY, 0.0, configuration->nBins_phiY,
                                       configuration->nBins_phiY, 0.0, configuration->nBins_phiY);

    calculateN1N1_H2H2H2(part1Histos->h_n1_phiY, part2Histos->h_n1_phiY, h_n1n1_phiYPhiY, 1.0, 1.0);
    calculateN1N1_H1H1H2(part1Histos->h_n1_pt, part2Histos->h_n1_pt, h_n1n1_ptPt, 1.0, 1.0);

    calculateN1N1_H2H2H2(part1Histos->h_spt_phiY, part2Histos->h_spt_phiY, h_pt1pt1_phiYPhiY, 1.0, 1.0);

    reduce(h_n1n1_phiYPhiY, h_n1n1_DyDphi, configuration->nBins_y, configuration->nBins_phi);
    reduce(h_pt1pt1_phiYPhiY, h_pt1pt1_DyDphi, configuration->nBins_y, configuration->nBins_phi);

    TF2* minusone = new TF2("minusone", "-1",
                            h_R2_DyDphi->GetXaxis()->GetBinLowEdge(1),
                            h_R2_DyDphi->GetXaxis()->GetBinUpEdge(h_R2_DyDphi->GetNbinsX()),
                            h_R2_DyDphi->GetYaxis()->GetBinLowEdge(1),
                            h_R2_DyDphi->GetYaxis()->GetBinUpEdge(h_R2_DyDphi->GetNbinsY()));
    double pt_1_average = part1Histos->h_pt_phiY->Integral() / (configuration->nBins_y * configuration->nBins_phi);
    double pt_2_average = part2Histos->h_pt_phiY->Integral() / (configuration->nBins_y * configuration->nBins_phi);
    if (reportInfo()) {
      cout << "calculateDerivedDiffHistograms(...) avgPt1 y = " << pt_1_average << endl;
      cout << "calculateDerivedDiffHistograms(...) avgPt2 Y = " << pt_2_average << endl;
    }

    /* calculate R2 */
    h_R2_DyDphi->Reset();
    h_R2_DyDphi->Add(h_n2_DyDphi, bincorrection);
    h_R2_DyDphi->Divide(h_n1n1_DyDphi);
    h_R2_DyDphi->Add(minusone);
    shiftY(*h_R2_DyDphi, *h_R2_DyDphi_shft, configuration->nBins_Dphi_shft);

    /* calculate G2 */
    h_G2_DyDphi->Reset();
    h_G2_DyDphi->Add(h_ptpt_DyDphi, bincorrection);
    h_G2_DyDphi->Add(h_pt1pt1_DyDphi, -1.0);
    h_G2_DyDphi->Divide(h_n1n1_DyDphi);
    h_G2_DyDphi->Scale(1.0 / (pt_1_average * pt_2_average));
    shiftY(*h_G2_DyDphi, *h_G2_DyDphi_shft, configuration->nBins_Dphi_shft);

    /* calculate P2 */
    h_P2_DyDphi->Reset();
    h_P2_DyDphi->Add(h_dptdpt_DyDphi, 1.0 / (pt_1_average * pt_2_average));
    h_P2_DyDphi->Divide(h_n2_DyDphi);
    shiftY(*h_P2_DyDphi, *h_P2_DyDphi_shft, configuration->nBins_Dphi_shft);

    /* calculate BF from R2 */
    double n1_2 = part1Histos->h_n1_phiY->Integral();
    double volume = kTWOPI * (part2Histos->h_n1_phiY->GetXaxis()->GetBinUpEdge(part2Histos->h_n1_phiY->GetNbinsX()) - part2Histos->h_n1_phiY->GetXaxis()->GetBinLowEdge(1));
    h_R2bf12_DyDphi_shft->Reset();
    h_R2bf12_DyDphi_shft->Add(h_R2_DyDphi_shft, n1_2 / volume);

    /* calculate Pratt's BF component from profiles */
    /* in principle the second single division should be automatically incorporated by the profiles when the same sign combination is substracted */
    /* we havo to compensate for the extra downscale by the profile mechanism and convert it to a density distribution */
    /* No, we don't! all the compensation should be taken care by the TProfile mechanism */
    p_PrattBf_DyDphi_shft->Reset();
    p_PrattBf_DyDphi_shft->Add(p_n2_DyDphi, 1.0);

    /* calculate Pratt's BF component from n2 and n1 histograms */
    TH2* tbf = (TH2*)h_n2_DyDphi->Clone(TString::Format("%s_clone", h_PrattBf_DyDphi_shft->GetName()).Data());
    tbf->Reset();
    tbf->Add(h_n2_DyDphi, h_n1n1_DyDphi, 1.0, -1.0);
    shiftY(*tbf, *h_PrattBf_DyDphi_shft, configuration->nBins_Dphi_shft);
    h_PrattBf_DyDphi_shft->Scale(1.0 / (n1_2 * h_PrattBf_DyDphi_shft->GetXaxis()->GetBinWidth(1) * h_PrattBf_DyDphi_shft->GetYaxis()->GetBinWidth(1)));
    delete tbf;

    delete h_n1n1_phiYPhiY;
    delete h_pt1pt1_phiYPhiY;
    delete minusone;
  }
  if (reportDebug())
    cout << "ParticlePairDerivedDiffHistos::calculateDerivedHistograms() completed." << endl;
}

ClassImp(ParticlePairDerivedDiffHistos)
