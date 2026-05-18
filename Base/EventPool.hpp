// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

#ifndef EventPool_hpp
#define EventPool_hpp

#include <TError.h>
#include <TMath.h>
#include "Particle.hpp"

using namespace std;

/////////////////////////////////////
// Class MiniParticle
/////////////////////////////////////
class MiniParticle
{
 public:
  MiniParticle()
    : pt(0),
      y(0),
      eta(0),
      phi(0),
      ixID(0),
      weight(1.0)
  {
  }

  ~MiniParticle()
  {
  }

  MiniParticle(const MiniParticle& other)
    : pt(other.pt),
      y(other.y),
      eta(other.eta),
      phi(other.phi),
      ixID(other.ixID),
      weight(other.weight)
  {
  }

  MiniParticle& operator=(const MiniParticle& other)
  {
    if (&other != this) {
      pt = other.pt;
      y = other.y;
      eta = other.eta;
      phi = other.phi;
      ixID = other.ixID;
      weight = other.weight;
    }
    return *this;
  }

  MiniParticle& operator=(const Particle& part)
  {
    pt = part.pt;
    y = part.y;
    eta = part.eta;
    phi = part.phi;
    ixID = part.ixID;
    weight = part.weight;

    return *this;
  }

  void printProperties(ostream& output)
  {
    output << "      pt: " << pt << endl;
    output << "       y: " << y << endl;
    output << "     eta: " << eta << endl;
    output << "     phi: " << phi << endl;
    output << "    ixID: " << ixID << endl;
    output << "  weight: " << weight << endl;
  }

  double getPx() { return pt * cos(phi); }
  double getPy() { return pt * sin(phi); }
  double getPz() { return pt / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-eta))); }
  double getE() { return 1.0 / TMath::TanH(y); }
  void getEPxPyPz(double* st)
  {
    st[0] = getE();
    st[1] = getPx();
    st[2] = getPy();
    st[3] = getPz();
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Data Members
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  double pt, y, eta;
  float phi;
  int ixID;      ///< the index of the filter which accepted the MiniParticle (kind of PID index)
  double weight; ///< per-particle correction weight (1.0 = no correction); inherited from Particle
};

/* Stores a pool of events in the shape of their particles       */
/* To avoid memory moves each event is a particle factory which  */
/* never gets deallocated, it is only reset to be fill with the  */
/* next event particles                                          */
/* The pool allocates an additional factory for the filling of   */
/* the new event at the same time the particles are selected     */
/* This additional factory is only part of the pool when the     */
/* user enqueue it                                               */
/* The pool should be full to start operating properly as a      */
/* mixed event pool. This means                                  */
/* - pool is full with the desired number of events but it       */
/* has the additional factory above described                    */
/* - pool acts as FIFO                                           */
/* - getNewEventStore gets the extra factory reseted             */
/*    - extra factory always at rear plus one                    */
/* - queue                                                       */
/*    - moves front to the next                                  */
/*    - moves rear to the next (actually the extra factory)      */
/* - get next with null returns the factory at front             */
/* - get next with position returns the next factory             */
/* - get next with rear factory returns null                     */
class EventPool
{
 public:
  EventPool(uint poolsize, uint factorysize)
    /* always an additional factory */
    : size(poolsize + 1),
      front(-1),
      rear(-1),
      pool{size, nullptr}
  {
    for (uint i = 0; i < size; ++i) {
      pool[i] = new Factory<MiniParticle>;
      pool[i]->initialize(factorysize);
    }
    Info("EventPool::EventPool()", "Built an event pool of %d events each with %d particles", size, factorysize);
  }

  ~EventPool()
  {
    for (auto ep : pool) {
      delete ep;
    }
    Info("EventPool::~EventPool()", "Destroyed an event pool of %d events", size);
  }

  bool isFull()
  {
    /* always an available factory after rear */
    //    Info("EventPool::isFull()", "Front: %d, rear: %d, so: %s", front, rear, ((rear + 2) % int(size) == front) ? "true" : "false");
    return ((rear + 2) % int(size) == front);
  }

  bool isEmpty()
  {
    //    Info("EventPool::isEmpty()", "Front: %d, rear: %d, so: %s", front, rear, (front == -1) ? "true" : "false");
    return (front == -1);
  }

  /* allways returns the available factory in rear plus one */
  Factory<MiniParticle>* getNewEventStore(int& ix)
  {
    ix = (rear + 1) % int(size);
    pool[ix]->reset();
    //    Info("EventPool::getNewEventStore()", "Front: %d, rear: %d, returning: %d", front, rear, ix);
    return pool[ix];
  }

  /* the user orders to enqueue the factory he got with the passed index */
  /* it is responsibility of the user to call the function once the      */
  /* factory has been properly filled                                    */
  /* The passed index has to be allways rear plus one                    */
  void enQueue(int ix)
  {
    //    Info("EventPool::enQueue()", "BEFORE: front: %d, rear: %d, enqueuing: %d", front, rear, ix);
    if (ix != (rear + 1) % int(size)) {
      Fatal("EventPool::enQueue", "Wrong passed index for enqueuing. Front: %d, rear: %d, ix: %d", front, rear, ix);
    }
    if (isFull()) {
      rear = (rear + 1) % int(size);
      front = (front + 1) % int(size);
    } else {
      if (front != -1 and front != 0) {
        Fatal("EventPool::enQueue", "Wrong behavior. Pool is not full with front: %d and rear: %d", front, rear);
      }
      front = 0;
      rear = (rear + 1) % int(size);
    }
    //    Info("EventPool::enQueue()", "AFTER: front: %d, rear: %d", front, rear);
  }

  Factory<MiniParticle>* getNextIndex(int& ix)
  {
    //    Info("EventPool::getNextIndex()", "BEFORE: front: %d, rear: %d, passing: %d", front, rear, ix);
    if (ix == -1) {
      ix = front;
      //      Info("EventPool::getNextIndex()", "AFTER: front: %d, rear: %d, returning: %d", front, rear, ix);
      return pool[front];
    } else {
      if (ix == rear) {
        ix = -1;
        //        Info("EventPool::getNextIndex()", "AFTER: front: %d, rear: %d, returning: %d and nullptr", front, rear, ix);
        return nullptr;
      } else {
        ix = (ix + 1) % int(size);
//        Info("EventPool::getNextIndex()", "AFTER: front: %d, rear: %d, returning: %d", front, rear, ix);
        return pool[ix];
      }
    }
  }

 private:
  uint size;
  int front;
  int rear;
  std::vector<Factory<MiniParticle>*> pool;
};

#endif /* EventPool_hpp */
