// -----------------------------------------------------------------------------
// Copyright (C) 2017  Ludovic LE FRIOUX
//
// This file is part of PaInleSS.
//
// PaInleSS is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
// -----------------------------------------------------------------------------

#pragma once

#include "../utils/Parameters.h"
#include "../utils/Threading.h"
#include "../working/SequentialWorker.h"
#include "../working/WorkingQueue.h"
#include "../working/WorkingStrategy.h"

#include <map>

using namespace std;

static void * mainMasterDivideAndConquer(void * arg);

class DivideAndConquer : public WorkingStrategy
{
public:
   DivideAndConquer();

   ~DivideAndConquer();

   void solve(const vector<int> & cube);

   void join(WorkingStrategy * strat, SatResult res,
             const vector<int> & model);
  
   void setInterrupt();

   void unsetInterrupt();

   void waitInterrupt();

   vector<int> getDivisionVariables(int k = 1);

   void setPhase(int var, bool value);

   void bumpVariableActivity(int var, int times);

protected:
   int cloneStrategy;
   int divisionStrategy;

   int numSplit;

   atomic<bool> strategyEnding;
   atomic<bool> waitJob;

   atomic<int> nDivisions; //number of divisions.
   atomic<int> nCancelledDivisions; //number of cancelled divisions.
   vector<double> timesLog;
   vector<double> splittingTimesLog;

   friend void * mainMasterDivideAndConquer(void * arg);

   Thread * master;

   pthread_mutex_t mutexStart;
   pthread_mutex_t mutexListsWorkers;
   pthread_mutex_t mutexFirstSolution;

   pthread_cond_t  mutexCondStart;
   pthread_cond_t mutexCondLists;

   atomic<int> nWorkers;

   vector<WorkingStrategy *> overs; //This vector contains WS that has finished working and ready to take another job
   vector<WorkingStrategy *> workers; //This vector contains WS that are working

   vector<WorkingStrategy *> workersSplitting; //This vector contains the workers that will have to split when they join.
   vector<WorkingStrategy *> oversSplitting;  //This vector contains some overs that will be split soon by a worker in workersSplitting.
     
   map<WorkingStrategy *,vector<int>> cubes; //This vector contains the cube for each worker
   map<WorkingStrategy *,double> times; //The vector contains the absolute time when each worker started solving its actual cube.

   vector<int> actualCube;

   WorkingQueue * workQueue;
};
