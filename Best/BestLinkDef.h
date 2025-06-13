#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class BestConfiguration + ;
#pragma link C++ class BestEventReader < AnalysisConfiguration::kRapidity> + ;
#pragma link C++ class BestEventReader < AnalysisConfiguration::kPseudorapidity> + ;
#pragma link C++ class BestAnalysisConfiguration + ;

#endif
