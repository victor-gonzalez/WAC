// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.                                   *
 * All rights reserved.                                                  *
 * Based on the ROOT package and environment                             *
 *                                                                       *
 * For the licensing terms see LICENSE.                                  *
 *************************************************************************/
/**
 \class GenericEventAnalyzer
 \ingroup WAC

 Utility class used to define graph parameters
 */

#include "TaskConfiguration.hpp"
ClassImp(TaskConfiguration);

TaskConfiguration::TaskConfiguration()
:
name("TaskConfiguration"),
type("Generic"),
version("1.0"),
loadHistograms(false),
createHistograms(true),
scaleHistograms(false),
calculateDerivedHistograms(false),
saveHistograms(true),
resetHistograms(false),
clearHistograms(false),
forceHistogramsRewrite(true),
inputPath("./"),
configurationFileName("configuration.txt"),
rootInputFileName(""),
outputPath("./"),
rootOuputFileName("output.root"),
dataInputPath(""),
dataInputFileName(""),
dataInputTreeName("tree"),
dataInputFileMinIndex(0),
dataInputFileMaxIndex(1)
{
/* no ops */
}

TaskConfiguration::TaskConfiguration(const TString & theName, const TString & theType, const TString & theVersion)
:
name(theName),
type(theType),
version(theVersion),
loadHistograms(false),
createHistograms(true),
scaleHistograms(false),
calculateDerivedHistograms(false),
saveHistograms(true),
resetHistograms(false),
clearHistograms(false),
forceHistogramsRewrite(true),
inputPath("./"),
configurationFileName("configuration.txt"),
rootInputFileName(""),
outputPath("./"),
rootOuputFileName("output.root")
{
/* no ops */
}

TaskConfiguration::~TaskConfiguration()
{

}

TString TaskConfiguration::getName() const
{
  return name;
}
TString TaskConfiguration::getType() const
{
  return type;
}
TString TaskConfiguration::getVersion() const
{
  return version;
}

void TaskConfiguration::printTaskConfiguration(ostream & os)
{
  os
  << " ---------------------------------------------" << endl
  << "               outputPath: " << name << endl
  << "    TaskConfigurationFileName: " << type << endl
  << "       rootOutputFileName: " << version << endl
  << " ---------------------------------------------" << endl;
}

TString TaskConfiguration::getInputRootFileName()
{
  TString inputName = inputPath;
  inputName += rootInputFileName;
  return inputName;
}

 TString TaskConfiguration::getOutputRootFileName()
{
  TString outputName = outputPath;
  outputName += rootOuputFileName;
  return outputName;
}
