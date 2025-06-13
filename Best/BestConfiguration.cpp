// Author: Claude Pruneau   09/25/2019

/*************************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 *************************************************************************/

#include "BestConfiguration.hpp"

ClassImp(BestConfiguration);

BestConfiguration::BestConfiguration()
  : TaskConfiguration("BEST", "EVREADER", "1.0")
{
}

////////////////////////////////////////////////////
// Print this configuration to the given stream
////////////////////////////////////////////////////
void BestConfiguration::printConfiguration(ostream& os)
{
  printTaskConfiguration(os);
  os
    << "    BEST   Parameters: " << endl
    << " ------------------------------------------------------------------------------------------" << endl;
}
