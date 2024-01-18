// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef WAC_Particle
#define WAC_Particle
#include "Factory.hpp"
using namespace std;

/////////////////////////////////////
// Class Particle
/////////////////////////////////////
class Particle
{
 public:
  Particle();
  ~Particle();
  Particle(const Particle& other);
  Particle& operator=(const Particle& other);
  void printProperties(ostream& output);
  void setPxPyPzE(double p_x, double p_y, double p_z, double p_e);

//  void setPidPxPyPzE(long pid,long charge, double p_x, double p_y, double p_z, double p_e);
//  void setPidPtPhiYEta(long _id, long _ch, double _pT, double _phi, double _y, double _eta);
  void setPidPxPyPzE(long pid, long charge, double p_x, double p_y, double p_z, double p_e, long pidMother=-1);
  void setPidPtPhiYEta(long _id, long _ch, double _pT, double _phi, double _y, double _eta, long pidMother=-1);


  void boost(double ax, double ay, double az);
  void boostRapidity(double boost);

  double getPx() { return px; }
  double getPy() { return py; }
  double getPz() { return pz; }
  double getE() { return e; }
  void getEPxPyPz(double* st)
  {
    st[0] = e;
    st[1] = px;
    st[2] = py;
    st[3] = pz;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  long pid;
  long charge;
  double px, py, pz, e, pt, y, eta;
  float phi;
  int ixYEtaPhi;
  int ixID; ///< the index of the filter which accepted the particle (kind of PID index)
  long pidMother;
  static int factorySize;
  static Factory<Particle>* factory;
  static Factory<Particle>* getFactory();
};

#endif /* WAC_Particle */
