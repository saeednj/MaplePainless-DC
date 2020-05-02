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

#include "../utils/Logger.h"
#include "../working/SequentialWorker.h"

#include <unistd.h>

using namespace std;

// Main executed by worker threads
void * mainWorker(void *arg)
{
   SequentialWorker * sq = (SequentialWorker *)arg;

   SatResult res = UNKNOWN;

   vector<int> model;

   while (globalEnding == false) {

      // Interrupt
      pthread_mutex_lock(&sq->interruptLock);

      if (sq->interrupt) {
         pthread_cond_wait(&sq->interruptCond, &sq->interruptLock);
      }

      pthread_mutex_unlock(&sq->interruptLock);


      // Solving phase
      sq->waitInterruptLock.lock();

      res = sq->solver->solve(sq->actualCube);

      sq->waitInterruptLock.unlock();

      // Report
      if (res == SAT) {
         model = sq->solver->getModel();
      }

      sq->join(NULL, res, model);

      model.clear();
   }

   return NULL;
}

// Constructor
SequentialWorker::SequentialWorker(SolverInterface * solver_)
{
   solver    = solver_;
   interrupt = true;

   pthread_mutex_init(&interruptLock, NULL);
   pthread_cond_init (&interruptCond, NULL);

   worker = new Thread(mainWorker, this);
}

// Destructor
SequentialWorker::~SequentialWorker()
{
   globalEnding = false;

   setInterrupt();
   waitInterrupt();

   worker->join();
   delete worker;

   pthread_mutex_destroy(&interruptLock);
   pthread_cond_destroy (&interruptCond);

   if(solver!=NULL)solver->release();
}

void
SequentialWorker::solve(const vector<int> & cube)
{
   actualCube = cube;
   
   unsetInterrupt();
}

void
SequentialWorker::join(WorkingStrategy * winner, SatResult res,
                       const vector<int> & model)
{
   if (globalEnding)
      return;

   if (parent == NULL) {
      globalEnding = true;
      finalResult  = res;

      if (res == SAT) {
         finalModel = model;
      }
   } else {
      parent->join(this, res, model);
   }
}

void
SequentialWorker::waitInterrupt()
{
   waitInterruptLock.lock();
   waitInterruptLock.unlock();
}

void
SequentialWorker::setInterrupt()
{
   interrupt = true;
   if(solver!=NULL)solver->setSolverInterrupt();
}

void
SequentialWorker::unsetInterrupt()
{
  if(solver!=NULL){solver->unsetSolverInterrupt();

   interrupt = false;

   pthread_mutex_lock  (&interruptLock);
   pthread_cond_signal (&interruptCond);
   pthread_mutex_unlock(&interruptLock);
  }
}

vector<int>
SequentialWorker::getDivisionVariables(int k)
{
   return solver->getDivisionVariables(k);
}


void
SequentialWorker::setPhase(const int var, const bool phase)
{
   solver->setPhase(var, phase);
}

void
SequentialWorker::bumpVariableActivity(const int var, const int times)
{
   solver->bumpVariableActivity(var, times);
}
