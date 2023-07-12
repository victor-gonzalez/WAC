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
#include <TBufferJSON.h>
#include "TRandom.h"
#include "../Base/AnalysisConfiguration.hpp"
#include "../Base/TwoPartDiffCorrelationAnalyzer.hpp"
#include "../Base/TwoPartDiffCorrelationAnalyzerME.hpp"
#include "../Pythia/PythiaAnalysisConfiguration.hpp"

const int ncomb = 2;
/* always charged particles first, positive and negative after it */
std::vector<std::string> partname;
const char* combname[ncomb] = {"CI", "CD"};
const int ncorrpart = 3;
/* balance function go its own way */
const char* corrfname[ncorrpart] = {"P2", "R2", "G2"};
/* balance function */
const int nbf = 8;
const char* bfnames[nbf] = {"R2BF", "R2BFPratt1bar2", "R2BFPrattbar12", "N2PrattBF", "N2PrattBF1bar2", "N2PrattBFbar12", "PrattBF1bar2", "PrattBFbar12"};

#define LONGITUDINAL AnalysisConfiguration::kRapidity

int nsamples = 10;

TFile* getSampleFile(PythiaAnalysisConfiguration* conf, int irap, int isample)
{
  TFile* f = nullptr;

  std::string filename = TString::Format("BUNCH%02d/Output/%s",
                                         isample + 1,
                                         TString::Format("%s.root", TString::Format(conf->outputfname.c_str(), int(conf->abs_y[irap] * 10)).Data()).Data())
                           .Data();
  f = new TFile(filename.c_str());
  if (f == nullptr or not f->IsOpen()) {
    Fatal("statUncertain::getSampleFile", "File %s not found. ABORTING!!!", filename.c_str());
  }
  return f;
}

TDirectory* getSampleDirectory(PythiaAnalysisConfiguration* conf, TFile* file, int irap, const char* cname)
{
  TDirectory* dir = nullptr;

  std::string dirname = TString::Format("%s%s_%s",
                                        TString::Format(conf->outputfname.c_str(), int(conf->abs_y[irap] * 10)).Data(),
                                        TString::Format(conf->taskname.c_str(), "Pairs").Data(),
                                        cname)
                          .Data();

  dir = (TDirectory*)file->GetDirectory(dirname.c_str());
  if (dir == nullptr) {
    Fatal("statUncertain::getSampleDirectory", "Directory %s not found in file %s. ABORTING!!!", dirname.c_str(), file->GetName());
  }
  return dir;
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

TList* extractSampleResults(Option_t* opt, PythiaAnalysisConfiguration* conf, TFile* samplefile, AnalysisConfiguration* ac, int irap, std::string evfltr, int isample)
{
  std::string feeddownrej = "none";
  char* cname = new char[64];
  auto getEvenFilter = [&cname](auto const& evfltr) {
    EventFilter* eventFilter = nullptr;
    if (evfltr == "MB") {
      sprintf(cname, "MB");
      eventFilter = new EventFilter(EventFilter::MinBias, 0, 0);
    } else {
      float min = 0;
      float max = 0;
      sscanf(evfltr.c_str(), "%f-%f", &min, &max);
      sprintf(cname, "C%dto%d", int(min), int(max));
      eventFilter = new EventFilter(EventFilter::Centrality, min, max);
    }
    return eventFilter;
  };

  /* we will control what goes to the directory tree */
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  std::vector<ParticleFilter<LONGITUDINAL>*> particleFilters;
  if (conf->inrapidity) {
    for (auto& part : conf->tpairs) {
      auto filter = PythiaAnalysisConfiguration::particleFilter<LONGITUDINAL>(part, feeddownrej, ac);
      if (filter != nullptr) {
        particleFilters.push_back(filter);
      } else {
        return nullptr;
      }
    }
  } else {
    Error("extractSampleResults", "Tasks in pseudorapidity still not supported by the analysis. Please, fix it!!");
    delete[] cname;
    return nullptr;
  }

  /* event selection at the analysis task level */
  EventFilter* eventFilter = getEvenFilter(evfltr);
  Event* event = Event::getEvent();

  /* the pairs taskname */
  TString taskName = TString::Format(conf->taskname.c_str(), TString::Format("PairsFDRej%s", feeddownrej.c_str()).Data());
  TwoPartDiffCorrelationAnalyzer<LONGITUDINAL>* eventanalyzer = new TwoPartDiffCorrelationAnalyzer<LONGITUDINAL>(taskName.Data(), ac, event, eventFilter, particleFilters);

  if (!TString(opt).Contains("verb"))
    eventanalyzer->setReportLevel(MessageLogger::Error);
  else
    eventanalyzer->setReportLevel(MessageLogger::Info);

  TDirectory* mydir = getSampleDirectory(conf, samplefile, irap, cname);
  if (mydir == nullptr)
    return nullptr;

  eventanalyzer->loadBaseHistograms(mydir);
  eventanalyzer->finalize();

  TList *list = new TList();
  list->SetOwner(kFALSE);

  /* the pair single histos */
  int npart = partname.size();
  bool doeta = TString(opt).Contains("eta");
  for (int ipart = 0; ipart < npart; ++ipart) {   /* first component of the pair */
    for (int jpart = 0; jpart < npart; ++jpart) { /* second component of the pair */
      TList* plist = new TList();                 /* a list per pair */
      plist->SetOwner(kTRUE);
      auto addToList = [plist, &cname, isample](auto h) {
        plist->Add(h->Clone(TString::Format("%s%s_Sub%02d", h->GetName(), cname, isample)));
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
      auto addToList = [plist, &cname, isample](auto h) {
        plist->Add(h->Clone(TString::Format("%s%s_Sub%02d", h->GetName(), cname, isample)));
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
      auto addToList = [plist, &cname, isample](auto h) {
        plist->Add(h->Clone(TString::Format("%s%s_Sub%02d", h->GetName(), cname, isample)));
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

  delete mydir;
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
  if (argc > 5 or argc < 5) {
    Fatal("main", "Wrong number of arguments. Use statUncertain prodtag option rapix evtfilterix");
  }

  std::string prodtag = argv[1];
  Option_t* opt = argv[2];
  int ixrap = stoi(argv[3]);
  int ixef = stoi(argv[4]);

  TTimeStamp now = TTimeStamp();
  if (!TString(opt).Contains("verb"))
    gErrorIgnoreLevel = kWarning;

  // =========================
  // Short configuration from json file
  // =========================
  PythiaAnalysisConfiguration* conf = nullptr;
  std::ifstream inf("configuration.json");
  if (inf.is_open()) {
    TString json;
    std::string line;
    while (getline(inf, line)) {
      json += line;
    }
    inf.close();
    /* let's produce the configuration object out of the string */
    TBufferJSON::FromJSON(conf, json);
    /* show it to chek it */
    TString test_json = TBufferJSON::ToJSON(conf);
    std::cout << test_json << std::endl;
  } else {
    printf("ERROR:Configuration file not found. ABORTING!!!\n");
    return 0;
  }

  /* store the selected particle species */
  partname = conf->tpairs;

  /* Cummulate errors by default */
  TH1::SetDefaultSumw2(kTRUE);

  auto getAnalysisConfiguration = [&conf, &ixrap]() {
    AnalysisConfiguration* ac = new AnalysisConfiguration("PYTHIA8", "UCM-miniWAC", "1.0");

    ac->outputPath = "./";
    ac->configurationFileName = "Config_DUKE.txt";
    ac->rootOuputFileName = TString::Format("Histograms%02dSub", nsamples).Data();

    // =========================
    // Short configuration
    // =========================
    float min_y = -conf->abs_y[ixrap];
    float max_y = conf->abs_y[ixrap];
    float min_pt = conf->min_pt;
    float max_pt = conf->max_pt;
    int nBins_y = int((max_y - min_y) / 0.1);

    ac->bin_edges_pt = conf->ptbins;
    ac->nBins_pt = conf->n_ptbins;
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

    return ac;
  };

  TFile* outputfile = new TFile(TString::Format("Histograms%02dSub_%s.root", nsamples, prodtag.c_str()), "UPDATE");
  char* ctitle = new char[64];
  auto getCentMultNames = [ctitle](auto const& evfltr) {
    if (evfltr == "MB") {
      sprintf(ctitle, "MB");
    } else {
      float min = 0;
      float max = 0;
      sscanf(evfltr.c_str(), "%f-%f", &min, &max);
      sprintf(ctitle, "%d-%d", int(min), int(max));
    }
  };

  /* the pair single histos, optionally the mixed events pair single histos, the balance function, plus the pair combined histos   */
  int npart = partname.size();
  int nmainlists = (TString(opt).Contains("me")) ? 2 * npart * npart + int(npart / 2) * int(npart / 2) + npart * (npart - 1) / 2 : npart * npart + int(npart / 2) * int(npart / 2) + npart * (npart - 1) / 2;
  std::vector<TObjArray> pairslists(nmainlists, TObjArray(nsamples));
  for (int ilst = 0; ilst < nmainlists; ++ilst) {
    pairslists[ilst].SetOwner(kTRUE);
  }

  for (Int_t isamp = 0; isamp < nsamples; isamp++) {
    TFile* samplefile = getSampleFile(conf, ixrap, isamp);
    if (samplefile != nullptr && samplefile->IsOpen()) {

      std::string ef = conf->teventfilter[ixef];
      Info("statUncertain", "Setup configuration for sample %d and event filter %s", isamp, ef.c_str());

      AnalysisConfiguration* ac = getAnalysisConfiguration();
      getCentMultNames(ef);

      Warning("statUncertain", "Processing sample %d for centrality %s", isamp, ctitle);
      TList* list = extractSampleResults(opt, conf, samplefile, ac, ixrap, ef, isamp);

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
    } else {
      Fatal("statUncertainPythia::main", "Cannot open rapidity %d sample file for sample %d. ABORTING!!!", ixrap, isamp);
    }
    samplefile->Close();
    delete samplefile;
  }

  /* we will control what goes to the directory tree */
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const char* raporeta = TString(opt).Contains("eta") ? "eta" : "y";
  outputfile->cd();
  for (int ipart = 0; ipart < npart; ++ipart) {
    for (int jpart = 0; jpart < npart; ++jpart) {
      int ilst = ipart * npart + jpart;
      TString pattern = TString::Format("Pythia8_%s%s%%s_D%sDphi_shft_%s", partname[ipart].c_str(), partname[jpart].c_str(), raporeta, ctitle);
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
        TString pattern = TString::Format("Pythia8_%s%s%%s_D%sDphi_shft_me_%s", partname[ipart].c_str(), partname[jpart].c_str(), raporeta, ctitle);
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
      TString pattern = TString::Format("Pythia8_%.2s%.2s%%s_D%sDphi_shft_%s", partname[ipart * 2].c_str(), partname[jpart * 2].c_str(), raporeta, ctitle);
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
      TString pattern = TString::Format("Pythia8_%s%s%%s_D%sDphi_shft_%s", partname[ipart].c_str(), partname[ipart + 1 + jpart].c_str(), raporeta, ctitle);
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

  outputfile->Close();
  delete[] ctitle;
  delete outputfile;
}
