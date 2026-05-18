// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/
/**
 \class Task
 \ingroup WAC

 Class defining Task
 */

#include "EventLoop.hpp"
ClassImp(EventLoop)

  EventLoop::EventLoop()
  /* bumped from 200 to 600 to accommodate the second analyzer set added when */
  /* detector effects are enabled (raw + reco), plus the generator and the    */
  /* DetectorEffectsTask itself                                               */
  : TaskCollection("EventLoop", nullptr, 600)
{
  if (reportDebug())
    cout << "EventLoop::CTOR Started" << endl;
}

EventLoop::~EventLoop()
{
  if (reportDebug())
    cout << "EventLoop::DTOR Started" << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Run the event loop as setup
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventLoop::run(long nEvents, int nReport)
{
  if (reportInfo("EventLoop", getName(), "run(...)"))
    cout << "Running for nEvents: " << nEvents << endl;
  postTaskOk();
  initialize();
  if (!isTaskOk()) {
    if (reportWarning("EventLoop", getName(), "run(...)"))
      cout << "Initialization failed. Abort." << endl;
    return;
  }
  for (long iEvent = 0; iEvent < nEvents; iEvent++) {
    if (isTaskOk())
      reset();
    if (isTaskOk())
      execute();
    if (isTaskOk() && iEvent % nReport == 0) {
      if (reportInfo("EventLoop", getName(), "run(...)"))
        cout << "Completed event # " << iEvent << endl;
    }
  }
  if (isTaskOk())
    finalize();

  if (reportInfo("EventLoop", getName(), "run(...)"))
    cout << "Completed nEvents: " << nEvents << endl;
  if (reportInfo("EventLoop", getName(), "run(...)"))
    cout << "Task finished with status : " << getTaskStatusName() << endl;
  if (reportDebug("EventLoop", getName(), "run(...)"))
    cout << "Completed." << endl;
}
