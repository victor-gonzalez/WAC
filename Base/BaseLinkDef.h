#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class HistogramCollection+;
#pragma link C++ class Histograms+;
#pragma link C++ class Task+;
#pragma link C++ class TaskCollection+;
#pragma link C++ class TaskConfiguration+;
#pragma link C++ class MessageLogger+;

#pragma link C++ class AnalysisConfiguration+;
#pragma link C++ class Event+;
#pragma link C++ class EventFilter+;
#pragma link C++ class EventHistos+;
#pragma link C++ class EventLoop+;
#pragma link C++ class GeneratorConfiguration+;
#pragma link C++ class TwoPartCorrelationAnalyzer < AnalysisConfiguration::kRapidity> + ;
#pragma link C++ class TwoPartCorrelationAnalyzer < AnalysisConfiguration::kPseudorapidity> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzer < AnalysisConfiguration::kRapidity, AnalysisConfiguration::kNoAdditionalOptions> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzer < AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillInvariantMass> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzer < AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillPratt> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzer < AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillPrattAndInvariantMass> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzer < AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kNoAdditionalOptions> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzer < AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillInvariantMass> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzer < AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillPratt> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzer < AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillPrattAndInvariantMass> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzerME < AnalysisConfiguration::kRapidity, AnalysisConfiguration::kNoAdditionalOptions> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzerME < AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillInvariantMass> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzerME < AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillPratt> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzerME < AnalysisConfiguration::kRapidity, AnalysisConfiguration::kFillPrattAndInvariantMass> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzerME < AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kNoAdditionalOptions> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzerME < AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillInvariantMass> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzerME < AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillPratt> + ;
#pragma link C++ class TwoPartDiffCorrelationAnalyzerME < AnalysisConfiguration::kPseudorapidity, AnalysisConfiguration::kFillPrattAndInvariantMass> + ;
#pragma link C++ class Factory < Particle> + ;
#pragma link C++ class Factory < MiniParticle> + ;
#pragma link C++ class Particle + ;
#pragma link C++ class MiniParticle + ;
#pragma link C++ class ParticleFilter < AnalysisConfiguration::kRapidity> + ;
#pragma link C++ class ParticleFilter < AnalysisConfiguration::kPseudorapidity> + ;
#pragma link C++ class ParticleHistos + ;
#pragma link C++ class ParticleAnalyzer < AnalysisConfiguration::kRapidity> + ;
#pragma link C++ class ParticleAnalyzer < AnalysisConfiguration::kPseudorapidity> + ;
#pragma link C++ class DetectorEffectsTask < AnalysisConfiguration::kRapidity> + ;
#pragma link C++ class DetectorEffectsTask < AnalysisConfiguration::kPseudorapidity> + ;
#pragma link C++ class ParticlePairCombinedHistos + ;
#pragma link C++ class ParticlePairBalanceFunctionDiffHistos + ;
#pragma link C++ class ParticlePairCombinedDiffHistos + ;
#pragma link C++ class ParticlePairDerivedHistos + ;
#pragma link C++ class ParticlePairDerivedDiffHistos + ;
#pragma link C++ class ParticlePairFilter + ;
#pragma link C++ class ParticlePairHistos+;
#pragma link C++ class ParticlePairDiffHistos + ;

#pragma link C++ class TrackAndPairConfiguration+;

#pragma link C++ class GeneralizedGaussianProfile+;
#pragma link C++ class GeneralizedGaussianConfiguration+;

#endif
