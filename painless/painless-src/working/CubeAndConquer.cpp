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

#include "../working/CubeAndConquer.h"
#include "../solvers/SolverInterface.h"
#include "../solvers/SolverFactory.h"
#include "../utils/Parameters.h"
#include "../working/SequentialWorker.h"
#include "../utils/Logger.h"
#include "../clauses/ClauseManager.h"
#include "../painless.h"

#include <unistd.h>
#include <algorithm>

using namespace std;

void * mainMasterCubeAndConquer(void * arg)
{
   CubeAndConquer * cc = (CubeAndConquer *)arg;

   vector<SequentialWorker *> targets;

   struct
   {
      bool operator()(SequentialWorker * a, SequentialWorker * b)
      {    
         return a->solver->getVariablesCount() < b->solver->getVariablesCount();
      }   
   } CompareSolverVariables;

   while(globalEnding == false) {
      // Wait until a job is notify by solve method
      pthread_mutex_lock(&cc->mutexStart);

      if (cc->waitJob == true) {
         pthread_cond_wait(&cc->mutexCondStart, &cc->mutexStart);
      }

      pthread_mutex_unlock(&cc->mutexStart);

      SolverInterface * cloneSolver =
         SolverFactory::cloneSolver(((SequentialWorker *)
                                     cc->slaves[0])->solver);

      if (sharers != NULL) {
         for (int i = 0; i < nSharers; i++) {
            sharers[i]->addConsumer(cloneSolver);
         }
      }

      SequentialWorker * cloneWorker = new SequentialWorker(cloneSolver);
      cc->workers.push_back(cloneWorker);
      cc->addSlave(cloneWorker);
     

      log(1, "Launch the portfolio thread\n");
      cc->slaves[0]->solve(cc->actualCube);
      
      log(1, "Launch the root\n");
      cloneWorker->solve(cc->actualCube);
      
      int time = 20;
      
      while (cc->strategyEnding == false) {
         log(1, "Next round in %d secondes\n", time);
         sleep(time);

         if (cc->strategyEnding)
            break;

         log(1, "Interrupt and wait workers\n");

         // Interrupt and wait workers
         for (size_t i=0; i<cc->workers.size(); i++) {
            cc->workers[i]->setInterrupt();
            cc->workers[i]->waitInterrupt();
         }

         log(1, "Erase %d UNSAT sub-problem(s)\n", cc->over.size());

         for (size_t i = 0; i < cc->over.size(); i++) {
            if (sharers != NULL) {
               for (int j = 0; j < nSharers; j++) {
                  sharers[j]->removeConsumer(cc->over[i]->solver);
               }
            }

            cc->workers.erase(remove(cc->workers.begin(), cc->workers.end(),
                     cc->over[i]), cc->workers.end());

            cc->slaves.erase(remove(cc->slaves.begin(), cc->slaves.end(),
                     cc->over[i]), cc->slaves.end());

            delete cc->over[i];
         }

         cc->over.clear();       

         targets.clear();

         if (cc->maxCpus  <= cc->workers.size()) {
            time *= 1.5;
         }

         int lowLimit = (int) cc->maxCpus / 4;

         if (cc->workers.size() <= lowLimit) {
            time = 20;
         }

         int actualWorkers = min(cc->maxCpus - 1, (int)cc->workers.size());

         for (size_t i = 0; i < actualWorkers; i++) {
            targets.push_back(cc->workers[i]);
         }

         log(1, "Divide %d sub-problems over %d\n", targets.size(),
               cc->workers.size());

         for (size_t i = 0; i < targets.size(); i++) {
            int var = targets[i]->getDivisionVariables()[0];

            log(1, "Target %d uses %d as division literal\n", i, var);

            SolverInterface * cloneSolver;
            cloneSolver = SolverFactory::cloneSolver(targets[i]->solver);

            ClauseExchange * cls = ClauseManager::allocClause(1);
            cls->lits[0]         = var;
            targets[i]->solver->addClause(cls);

            ClauseExchange * clsClone = ClauseManager::allocClause(1);
            clsClone->lits[0]         = -var;
            cloneSolver->addClause(clsClone);

            if (sharers != NULL) {
               for (int j=0; j < nSharers; j++) {
                  sharers[j]->addConsumer(cloneSolver);
               }
            }

            SequentialWorker * cloneWorker = new SequentialWorker(cloneSolver);
            cc->workers.push_back(cloneWorker);
            cc->addSlave(cloneWorker);

            cc->nJobs++;
         }

         if (cc->workers.size() * 2 > cc->maxNodes) {
            sort(cc->workers.begin(), cc->workers.end(),
                  CompareSolverVariables);
         }

         actualWorkers = min(cc->maxCpus - 1, (int)cc->workers.size());

         log(1, "Launch %d workers for this round\n", actualWorkers);

         for (size_t i = 0; i < actualWorkers; i++) {
            cc->workers[i]->solve(cc->actualCube);
         }

         cc->next = actualWorkers;
      }         
      
      cc->waitJob = true;
   } 

   return NULL;
}

CubeAndConquer::CubeAndConquer(int maxCpus_)
{
   maxCpus = maxCpus_;
   maxNodes = 8 * maxCpus;

   pthread_mutex_init(&mutexStart, NULL);
   pthread_cond_init (&mutexCondStart, NULL);

   waitJob = true;

   master = new Thread(mainMasterCubeAndConquer, this);
}

CubeAndConquer::~CubeAndConquer()
{
   master->join();
   delete master;

   pthread_mutex_destroy(&mutexStart);
   pthread_cond_destroy (&mutexCondStart);
}

   void
CubeAndConquer::solve(const vector<int> & cube)
{
   strategyEnding = false;
   nJobs          = 1;

   waitJob = false;

   pthread_mutex_lock (&mutexStart);
   pthread_cond_signal (&mutexCondStart);
   pthread_mutex_unlock(&mutexStart);
}

   void
CubeAndConquer::join(WorkingStrategy * strat, SatResult res,
      const vector<int> & model)
{
   if (res == UNKNOWN || strategyEnding)
      return;

   if (res == UNSAT && strat != slaves[0]) {
      nJobs--;

      if (nJobs.load() > 0) {
         log(1 , "UNSAT sub-problem resolved, %d left\n", nJobs.load());

         over.push_back((SequentialWorker *)strat);

         if (next.load() < workers.size()) {
            workers[next]->solve(actualCube);
            next++;
         }

         return;
      }
   }

   strategyEnding = true;

   setInterrupt();

   if (parent == NULL) { // If it is the top strategy
      globalEnding = true;
      finalResult  = res;

      SequentialWorker * winner = (SequentialWorker *)strat;
      log(0, "Winner: %d\n", winner->solver->id);

      if (res == SAT) {
         finalModel = model;
      }
   } else { // Else forward the information to the parent strategy
      parent->join(this, res, model);  
   }
}

   void
CubeAndConquer::setInterrupt()
{
   for (size_t i = 0; i < slaves.size(); i++) {
      slaves[i]->setInterrupt();
   }
}

   void
CubeAndConquer::unsetInterrupt()
{
   for (size_t i = 0; i < slaves.size(); i++) {
      slaves[i]->unsetInterrupt();
   }
}

   void
CubeAndConquer::waitInterrupt()
{
   for (size_t i = 0; i < slaves.size(); i++) {
      slaves[i]->waitInterrupt();
   }
}

   vector<int>
CubeAndConquer::getDivisionVariables(int k)
{
   return vector<int>(1, 0);
}

   void
CubeAndConquer::setPhase(const int var, const bool phase)
{
}

   void
CubeAndConquer::bumpVariableActivity(const int var, const int times)
{
}
