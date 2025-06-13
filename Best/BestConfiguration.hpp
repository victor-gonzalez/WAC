// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
#ifndef BESTCONFIGURATION
#define BESTCONFIGURATION
#include "TaskConfiguration.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////
// BEST Event reader configuration
////////////////////////////////////////////////////////////////////////////////////////////////////////
class BestConfiguration : public TaskConfiguration
{
 public:
  BestConfiguration();
  virtual ~BestConfiguration() {}
  void printConfiguration(ostream& os);

  ////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////

  ClassDef(BestConfiguration, 0)
};

#endif // BESTCONFIGURATION
