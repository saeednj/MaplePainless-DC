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

#include "../utils/Threading.h"
#include "../working/SequentialWorker.h"
#include "../working/WorkingStrategy.h"

using namespace std;

static void * mainMasterCubeAndConquer(void * arg);

class CubeAndConquer : public WorkingStrategy
{
public:
   CubeAndConquer(int maxCpus);

   ~CubeAndConquer();

   void solve(const vector<int> & cube);

   void join(WorkingStrategy * strat, SatResult res, const vector<int> & model);

   void setInterrupt();

   void unsetInterrupt();

   void waitInterrupt();

   vector<int> getDivisionVariables(int k = 1);

   void setPhase(const int var, const bool phase);

   void bumpVariableActivity(const int var, const int times);

protected:
   friend void * mainMasterCubeAndConquer(void * arg); 

   atomic<bool> strategyEnding;
   
   atomic<bool> waitJob;

   vector<int> actualCube;

   Thread * master;
   
   int maxNodes;
   
   int maxCpus;

   pthread_mutex_t mutexStart;
   pthread_cond_t  mutexCondStart;
   
   atomic<int> nJobs;
   atomic<int> next;

   vector<SequentialWorker *> over;
   vector<SequentialWorker *> workers;
};
