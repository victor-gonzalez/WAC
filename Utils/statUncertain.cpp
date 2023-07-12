//
//  statUncertain.C
//  CodeV3
//
//  Created by Claude Pruneau on 8/22/17.
//  Copyright © 2017 Claude Pruneau. All rights reserved.
//

#include <stdio.h>
#include <cstring>

//
//  cm_figure9_6.c
//  MyMC
//
//  Created by Claude Pruneau on 4/10/17.
//  Copyright © 2017 Claude Pruneau. All rights reserved.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include <TError.h>
#include <TObjArray.h>
#include <TTimeStamp.h>
#include <TList.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMath.h>
#include "TRandom.h"
#include "../Base/AnalysisConfiguration.hpp"
#include "../Base/TwoPartDiffCorrelationAnalyzer.hpp"
#include "../Base/TwoPartDiffCorrelationAnalyzerME.hpp"

const int npart = 6;
const int ncomb = 2;
/* always charged particles first, positive and negative after it */
const char* partname[npart] = {"PiP", "PiM", "KaP", "KaM", "PrP", "PrM"};
const char* combname[ncomb] = {"CI", "CD"};
const int ncorrpart = 3;
/* balance function go its own way */
const char* corrfname[ncorrpart] = {"P2", "R2", "G2"};
/* balance function */
const int nbf = 8;
const char* bfnames[nbf] = {"R2BF", "R2BFPratt1bar2", "R2BFPrattbar12", "N2PrattBF", "N2PrattBF1bar2", "N2PrattBFbar12", "PrattBF1bar2", "PrattBFbar12"};

#define LONGITUDINAL AnalysisConfiguration::kRapidity

std::vector<std::string> centfname = {
  "MB"};
int ncent = centfname.size();

std::vector<std::string> samplesfnames = {
  "BUNCH01/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH02/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH03/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH04/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH05/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH06/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH07/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH08/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH09/Output/PYTHIA8_PiKaPr_%s.root",
  "BUNCH10/Output/PYTHIA8_PiKaPr_%s.root",
};
int nsamples = samplesfnames.size();

TFile *getSampleFile(int icent,int isample) {
  TFile *f = nullptr;

  std::string filename = TString::Format(samplesfnames[isample].c_str(), centfname[icent].c_str()).Data();
  f = new TFile(filename.c_str());
  if (f == nullptr or not f->IsOpen()) {
    Error("statUncertain::getSampleFile","File %s not found. ABORTING!!!", filename.c_str());
  }
  return f;
}


TH2 *extractHistoMeanAndStDevFromSubSets(const TObjArray &listsarray, Int_t ih, const TString &name) {

  TH2 *h2 = dynamic_cast<TH2*>(((TList*) listsarray[0])->At(ih));

  if (h2 != NULL) {
    /* extract the name of the result histogram */
    TString title = h2->GetTitle();
    TString xtitle = h2->GetXaxis()->GetTitle();
    TString ytitle = h2->GetYaxis()->GetTitle();
    TString ztitle = h2->GetZaxis()->GetTitle();

    /* let's first create the profile we use as support */
    TProfile2D *p = new TProfile2D("auxprofile","auxprofile",
        h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(),
        h2->GetNbinsY(), h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax());

    for (Int_t ix = 0; ix < h2->GetNbinsX(); ix++) {
      Double_t x = h2->GetXaxis()->GetBinCenter(ix+1);
      for (Int_t iy = 0; iy < h2->GetNbinsY(); iy++) {
        Double_t y = h2->GetYaxis()->GetBinCenter(iy+1);
        for (Int_t iset = 0; iset < listsarray.GetEntries(); iset++) {
          p->Fill(x,y,((TH2*)(((TList*) listsarray[iset])->At(ih)))->GetBinContent(ix+1,iy+1));
        }
      }
    }
    h2 = p->ProjectionXY(name.Data());
    h2->SetTitle(title.Data());
    h2->SetXTitle(xtitle.Data());
    h2->SetYTitle(ytitle.Data());
    h2->SetZTitle(ztitle.Data());

    /* set the proper number of entries */
    Double_t nentries = 0;
    for (Int_t iset = 0; iset < listsarray.GetEntries(); iset++) {
      nentries += ((TH2*)(((TList*) listsarray[iset])->At(ih)))->GetEntries();
    }
    h2->SetEntries(nentries);

    delete p;
    return h2;
  }
  else {
    Error("extractHistoMeanAndStDevFromSubSets", "Wrong histogram type");
    return NULL;
  }
}

TList* extractMeanAndStDevFromSubSets(const TObjArray& listsarray, const TString& pattern, const char* name[])
{
  /* basically we receive an array of histograms lists */
  /* each array item corresponds to a results subset   */
  /* we extract an equivalent list of histograms which */
  /* corresponds to the mean of the subsets histograms */
  /* and have as errors the standard deviation from    */
  /* the mean on each bin                              */

  /* first, some consistency checks */
  Int_t nhistos = ((TList*) listsarray[0])->GetEntries();
  if (nhistos != ncorrpart and nhistos != (ncorrpart * ncomb) and nhistos != nbf)
    Error("extractMeanAndStDevFromSubSets", "Inconsistent number of histograms to average");
  for (Int_t iset = 0; iset < listsarray.GetEntries(); iset++) {
    if (nhistos != ((TList*) listsarray[iset])->GetEntries()) {
      Error("extractMeanAndStDevFromSubSets", "Inconsistent array of lists");
      return NULL;
    }
  }

  TList *list = new TList();
  list->SetOwner(kTRUE);
  for (Int_t ih = 0; ih < nhistos; ih++) {
    list->Add(extractHistoMeanAndStDevFromSubSets(listsarray,ih,TString::Format(pattern.Data(),name[ih])));
  }

  return list;
}

TList* extractSampleResults(Option_t* opt, AnalysisConfiguration* ac, int icent, int isample)
{

  /* we will control what goes to the directory tree */
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  std::vector<ParticleFilter<LONGITUDINAL>*> particleFilters;
  particleFilters.push_back(new ParticleFilter<LONGITUDINAL>(ParticleFilter<LONGITUDINAL>::Pion, ParticleFilter<LONGITUDINAL>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<LONGITUDINAL>(ParticleFilter<LONGITUDINAL>::Pion, ParticleFilter<LONGITUDINAL>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<LONGITUDINAL>(ParticleFilter<LONGITUDINAL>::Kaon, ParticleFilter<LONGITUDINAL>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<LONGITUDINAL>(ParticleFilter<LONGITUDINAL>::Kaon, ParticleFilter<LONGITUDINAL>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<LONGITUDINAL>(ParticleFilter<LONGITUDINAL>::Proton, ParticleFilter<LONGITUDINAL>::Positive, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));
  particleFilters.push_back(new ParticleFilter<LONGITUDINAL>(ParticleFilter<LONGITUDINAL>::Proton, ParticleFilter<LONGITUDINAL>::Negative, ParticleFilter<AnalysisConfiguration::kRapidity>::None, ac->min_pt, ac->max_pt, ac->min_y, ac->max_y));

  EventFilter* eventFilter = new EventFilter(EventFilter::MinBias, 0, 0);

  Event *event = Event::getEvent();

  TwoPartDiffCorrelationAnalyzer<LONGITUDINAL>* eventanalyzer = new TwoPartDiffCorrelationAnalyzer<LONGITUDINAL>("NarrowPiKaPrLa", ac, event, eventFilter, particleFilters);

  if (!TString(opt).Contains("verb"))
    eventanalyzer->setReportLevel(MessageLogger::Error);
  else
    eventanalyzer->setReportLevel(MessageLogger::Info);

  TFile* myfile = getSampleFile(icent, isample);
  if (myfile == nullptr) return nullptr;

  eventanalyzer->loadBaseHistograms(myfile);
  eventanalyzer->finalize();

  TList *list = new TList();
  list->SetOwner(kFALSE);

  /* the pair single histos */
  bool doeta = TString(opt).Contains("eta");
  for (int ipart = 0; ipart < npart; ++ipart) {   /* first component of the pair */
    for (int jpart = 0; jpart < npart; ++jpart) { /* second component of the pair */
      TList* plist = new TList();                 /* a list per pair */
      plist->SetOwner(kTRUE);
      auto addToList = [plist, icent, isample](auto h) {
        plist->Add(h->Clone(TString::Format("%s%s_Sub%02d", h->GetName(), centfname[icent].c_str(), isample)));
      };
      if (doeta) {
        addToList(eventanalyzer->pairs_Histos[ipart][jpart]->h_P2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_Histos[ipart][jpart]->h_R2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_Histos[ipart][jpart]->h_G2_DetaDphi_shft);
      } else {
        addToList(eventanalyzer->pairs_Histos[ipart][jpart]->h_P2_DyDphi_shft);
        addToList(eventanalyzer->pairs_Histos[ipart][jpart]->h_R2_DyDphi_shft);
        addToList(eventanalyzer->pairs_Histos[ipart][jpart]->h_G2_DyDphi_shft);
      }
      list->Add(plist);
    }
  }

#ifdef MIXED_EVENTS /* TODO: we have to template this */
  if (TString(opt).Contains("me")) {
    /* the pair single histos from ME */
    for (int ipart = 0; ipart < npart; ++ipart) {   /* first component of the pair */
      for (int jpart = 0; jpart < npart; ++jpart) { /* second component of the pair */
        TList* plist = new TList();                 /* a list per pair */
        plist->SetOwner(kTRUE);
        for (int icf = 0; icf < ncorrpart; ++icf) { /* the individual pair cf */
          TH2* h2 = NULL;
          switch (icf) {
            case 0: /* P2 */
              h2 = eventanalyzer->pairs_Histos_me[ipart][jpart]->h_P2_DetaDphi_shft;
              break;
            case 1: /* R2 */
              h2 = eventanalyzer->pairs_Histos_me[ipart][jpart]->h_R2_DetaDphi_shft;
              break;
            case 2: /* G2 */
              h2 = eventanalyzer->pairs_Histos_me[ipart][jpart]->h_G2_DetaDphi_shft;
              break;
            default:
              Fatal("extractSampleResults", "Wrong correlator index");
          }
          plist->Add(h2->Clone(TString::Format("%s%s_Sub%02d", h2->GetName(), centfname[icent].c_str(), isample)));
        }
        list->Add(plist);
      }
    }
  }
#endif

  /* the pair balance functions */
  for (int ipart = 0; ipart < int(npart / 2); ++ipart) {
    for (int jpart = 0; jpart < int(npart / 2); ++jpart) {
      TList* plist = new TList(); /* a list per pair */
      plist->SetOwner(kTRUE);
      auto addToList = [plist, icent, isample](auto h) {
        plist->Add(h->Clone(TString::Format("%s%s_Sub%02d", h->GetName(), centfname[icent].c_str(), isample)));
      };
      if (doeta) {
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_R2BF_DetaDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_R2BF_Pratt_1bar2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_R2BF_Pratt_bar12_DetaDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_PrattBF_DetaDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_PrattBF_1bar2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_PrattBF_bar12_DetaDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->p_PrattBF_1bar2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->p_PrattBF_bar12_DetaDphi_shft);
      } else {
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_R2BF_DyDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_R2BF_Pratt_1bar2_DyDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_R2BF_Pratt_bar12_DyDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_PrattBF_DyDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_PrattBF_1bar2_DyDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->h_PrattBF_bar12_DyDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->p_PrattBF_1bar2_DyDphi_shft);
        addToList(eventanalyzer->pairs_BFHistos[ipart][jpart]->p_PrattBF_bar12_DyDphi_shft);
      }
      list->Add(plist);
    }
  }

  /* the pair combined histos */
  for (int ipart = 0; ipart < npart; ++ipart) {                 /* first component of the pair */
    for (int jpart = 0; jpart < npart - (ipart + 1); ++jpart) { /* second component of the pair */
      TList* plist = new TList();                               /* a list per pair */
      plist->SetOwner(kTRUE);
      auto addToList = [plist, icent, isample](auto h) {
        plist->Add(h->Clone(TString::Format("%s%s_Sub%02d", h->GetName(), centfname[icent].c_str(), isample)));
      };
      if (doeta) {
        addToList(eventanalyzer->pairs_CIHistos[ipart][jpart]->h_P2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_CDHistos[ipart][jpart]->h_P2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_CIHistos[ipart][jpart]->h_R2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_CDHistos[ipart][jpart]->h_R2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_CIHistos[ipart][jpart]->h_G2_DetaDphi_shft);
        addToList(eventanalyzer->pairs_CDHistos[ipart][jpart]->h_G2_DetaDphi_shft);
      } else {
        addToList(eventanalyzer->pairs_CIHistos[ipart][jpart]->h_P2_DyDphi_shft);
        addToList(eventanalyzer->pairs_CDHistos[ipart][jpart]->h_P2_DyDphi_shft);
        addToList(eventanalyzer->pairs_CIHistos[ipart][jpart]->h_R2_DyDphi_shft);
        addToList(eventanalyzer->pairs_CDHistos[ipart][jpart]->h_R2_DyDphi_shft);
        addToList(eventanalyzer->pairs_CIHistos[ipart][jpart]->h_G2_DyDphi_shft);
        addToList(eventanalyzer->pairs_CDHistos[ipart][jpart]->h_G2_DyDphi_shft);
      }
      list->Add(plist);
    }
  }

  myfile->Close();

  delete myfile;
  delete eventanalyzer;
  delete eventFilter;
  for (auto pf : particleFilters) {
    delete pf;
  }

  /* back to normal */
  TH1::AddDirectory(oldstatus);

  return list;
}

/////////////////////////////////////////////////////////////////////////////////////////
// produce results with statistic uncertainties out of a set of sub-samples
/////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  if (argc > 2 or argc < 1) {
    Fatal("main", "Wrong number of arguments. Use statUncertain option");
  }

  Option_t* opt = argv[1];

  TTimeStamp now = TTimeStamp();
  if (!TString(opt).Contains("verb"))
    gErrorIgnoreLevel = kWarning;

  /* Cummulate errors by default */
  TH1::SetDefaultSumw2(kTRUE);

  TFile* outputfile = new TFile(TString::Format("Histograms%02dSub_%d_%d.root", nsamples, now.GetDate(), now.GetTime()), "RECREATE");

  for (int icent = 0; icent < ncent; icent++) {

    Info("statUncertain", "Setup configuration");

    AnalysisConfiguration* ac = new AnalysisConfiguration("PYTHIA8", "UCM-miniWAC", "1.0");

    ac->outputPath = "./";
    ac->configurationFileName = "Config_DUKE.txt";
    ac->rootOuputFileName = TString::Format("Histograms%02dSub", nsamples).Data();

    cout << "================================================================" << endl;
    cout << "" << endl;
    cout << "         running for = " << ac->getName() << endl;
    cout << "" << endl;
    cout << "================================================================" << endl;

    // =========================
    // Short configuration
    // =========================
    float min_y = -4;
    float max_y = 4;
    float min_pt = 0.0;
    float max_pt = 1e6;
    int nBins_y = int((max_y - min_y) / 0.1);
    int nBins_pt = int((max_pt - min_pt) / 0.1);

    ac->bin_edges_pt = {0.0,
                        0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                        0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4,
                        1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                        3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 13.0, 20.0};
    ac->nBins_pt = ac->bin_edges_pt.size() - 1;
    ac->min_pt = min_pt;
    ac->max_pt = max_pt;
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
    ac->fillPratt = false;
    ac->fillYorEta = LONGITUDINAL;

    ac->scaleHistograms = true;
    ac->createHistograms = false;
    ac->loadHistograms = true;
    ac->saveHistograms = false;
    ac->forceHistogramsRewrite = false;
    ac->calculateDerivedHistograms = true;

    /* the pair single histos, optionally the mixed events pair single histos, the balance function, plus the pair combined histos   */
    int nmainlists = (TString(opt).Contains("me")) ? 2 * npart * npart + int(npart / 2) * int(npart / 2) + npart * (npart - 1) / 2 : npart * npart + int(npart / 2) * int(npart / 2) + npart * (npart - 1) / 2;
    std::vector<TObjArray> pairslists(nmainlists, TObjArray(nsamples));
    for (int ilst = 0; ilst < nmainlists; ++ilst) {
      pairslists[ilst].SetOwner(kTRUE);
    }

    for (Int_t isamp = 0; isamp < nsamples; isamp++) {
      Warning("statUncertain", "Processing sample %d for centrality %s", isamp, centfname[icent].c_str());
      TList* list = extractSampleResults(opt, ac, icent, isamp);

      for (Int_t ilst = 0; ilst < nmainlists; ++ilst) {
        pairslists[ilst][isamp] = list->At(ilst);
      }
      /* we write the individual results if required */
      if (TString(opt).Contains("savesub")) {
        outputfile->cd();
        for (Int_t ixh = 0; ixh < list->GetEntries(); ixh++) {
          list->At(ixh)->Write();
        }
      }
      delete list;
    }

    /* we will control what goes to the directory tree */
    Bool_t oldstatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    const char* raporeta = TString(opt).Contains("eta") ? "eta" : "y";
    outputfile->cd();
    for (int ipart = 0; ipart < npart; ++ipart) {
      for (int jpart = 0; jpart < npart; ++jpart) {
        int ilst = ipart * npart + jpart;
        TString pattern = TString::Format("Pythia8_%s%s%%s_D%sDphi_shft_%s", partname[ipart], partname[jpart], raporeta, centfname[icent].c_str());
        TList* meanhlist = extractMeanAndStDevFromSubSets(pairslists[ilst], pattern, corrfname);
        for (int ixh = 0; ixh < meanhlist->GetEntries(); ixh++) {
          meanhlist->At(ixh)->Write();
        }
        delete meanhlist;
      }
    }
    if (TString(opt).Contains("me")) {
      int ilst = npart * npart;
      for (int ipart = 0; ipart < npart; ++ipart) {
        for (int jpart = 0; jpart < npart; ++jpart) {
          TString pattern = TString::Format("Pythia8_%s%s%%s_D%sDphi_shft_me_%s", partname[ipart], partname[jpart], raporeta, centfname[icent].c_str());
          TList* meanhlist = extractMeanAndStDevFromSubSets(pairslists[ilst++], pattern, corrfname);
          for (int ixh = 0; ixh < meanhlist->GetEntries(); ixh++) {
            meanhlist->At(ixh)->Write();
          }
          delete meanhlist;
        }
      }
    }

    /* the charge balance function */
    int ilst = npart * npart;
    if (TString(opt).Contains("me")) {
      ilst *= 2;
    }
    for (int ipart = 0; ipart < int(npart / 2); ++ipart) {
      for (int jpart = 0; jpart < int(npart / 2); ++jpart) {
        TString pattern = TString::Format("Pythia8_%.2s%.2s%%s_D%sDphi_shft_%s", partname[ipart * 2], partname[jpart * 2], raporeta, centfname[icent].c_str());
        TList* meanhlist = extractMeanAndStDevFromSubSets(pairslists[ilst++], pattern, bfnames);
        for (int ixh = 0; ixh < meanhlist->GetEntries(); ixh++) {
          meanhlist->At(ixh)->Write();
        }
        delete meanhlist;
      }
    }

    /* the CI, CD, combinations */
    const char* cfnamecomb[ncorrpart * ncomb] = {nullptr};
    for (int i = 0; i < ncorrpart; ++i) {
      for (int j = 0; j < ncomb; ++j) {
        char* buffer;
        cfnamecomb[i * ncomb + j] = buffer = new char[std::strlen(corrfname[i]) + strlen(combname[j]) + 1];
        sprintf(buffer, "%s%s", corrfname[i], combname[j]);
      }
    }
    /* we keep tracking with the previous ilst content */
    for (int ipart = 0; ipart < npart; ++ipart) {
      for (int jpart = 0; jpart < npart - (ipart + 1); ++jpart) {
        TString pattern = TString::Format("Pythia8_%s%s%%s_D%sDphi_shft_%s", partname[ipart], partname[ipart + 1 + jpart], raporeta, centfname[icent].c_str());
        TList* meanhlist = extractMeanAndStDevFromSubSets(pairslists[ilst++], pattern, cfnamecomb);
        for (int ixh = 0; ixh < meanhlist->GetEntries(); ixh++) {
          meanhlist->At(ixh)->Write();
        }
        delete meanhlist;
      }
    }

    for (auto aname : cfnamecomb) {
      delete aname;
    }
    /* back to normal */
    TH1::AddDirectory(oldstatus);
  }
  outputfile->Close();
  delete outputfile;
}
