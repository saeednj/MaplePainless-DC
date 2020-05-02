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

#include "painless.h"

#include "clauses/ClauseManager.h"
#include "sharing/SimpleSharing.h"
#include "sharing/HordeSatSharing.h"
#include "sharing/Sharer.h"
#include "solvers/SolverFactory.h"
#include "utils/Logger.h"
#include "utils/Parameters.h"
#include "utils/SatUtils.h"
#include "utils/System.h"
#include "working/SequentialWorker.h"
#include "working/Portfolio.h"
#include "working/CubeAndConquer.h"
#include "working/DivideAndConquer.h"

#include <unistd.h>

using namespace std;


// -------------------------------------------
// Declaration of global variables
// -------------------------------------------
atomic<bool> globalEnding(false);

Sharer ** sharers = NULL;

int nSharers = 0;

WorkingStrategy * working = NULL;

SatResult finalResult = UNKNOWN;

vector<int> finalModel;


// -------------------------------------------
// Main du framework
// -------------------------------------------
int main(int argc, char ** argv)
{
   Parameters::init(argc, argv);
   
   if (Parameters::getFilename() == NULL ||
       Parameters::isSet("h"))
   {
      printf("USAGE: %s [parameters] input.cnf\n", argv[0]);
      printf("Parameters:\n");
      printf("\t-solver=type\t\t glucose, lingeling, minisat, maple, " \
             "maplesat, or combo\n");
      printf("\t-lbd-limit=<INT>\t lbd limit of exported clauses\n");
      printf("\t-d=0...7\t\t diversification 0=none, 1=sparse, 2=dense," \
             " 3=random, 4=native, 5=1&4, 6=sparse-random, 7=6&4," \
             " default is 0.\n");
      printf("\t-c=<INT>\t\t number of cpus, default is 4.\n");
      printf("\t-wkr-strat=1...4\t 1=portfolio, 2=cube and conquer," \
             " 4=divide and conquer, default is portfolio\n");
      printf("\t-shr-strat=1...2\t 1=alltoall, 2=hordesat sharing," \
             " default is 0\n");
      printf("\t-shr-sleep=<INT>\t time in usecond a sharer sleep each" \
             " round, default 500000 (0.5s)\n");
      printf("\t-shr-lit=<INT>\t\t number of literals shared per round by the" \
             " hordesat strategy, default is 1500\n");
      printf("\t-no-model\t\t won't print the model if the problem is SAT\n");
      printf("\t-t=<INT>\t\t timeout in second, default is no limit\n");
      printf("\t-split-heur=1...3\t for D&C: splitting heuristic," \
	          " 1=VSIDS, 2=#flips, 3=propagation rate, default is 1\n");
      printf("\t-copy-mode=1...2\t for D&C: method to allocate new subspaces " \
             "to workers, 1=no copy: reuses the old solver, 2=clone: clones " \
             "the solver and delete old solver, default is 1\n");
      printf("\t-num-split=<INT>\t for D&C: number of splitting variables picked at" \
             "each splitting point, default is 1\n");

      return 0;
   }

   Parameters::printParams();

   setVerbosityLevel(Parameters::getIntParam("v", 0));

   int cpus = Parameters::getIntParam("c", 4);

   srand(time(NULL));

   // Create solvers
   vector<SolverInterface *> solvers;
   
   const string solverType = Parameters::getParam("solver");
   const int wkr_strat=Parameters::getIntParam("wkr-strat", 1);

   int nSolvers = cpus;
   if (wkr_strat == 2 || wkr_strat == 5 || (wkr_strat == 4 &&
            Parameters::getIntParam("copy-mode", 1) == 2)){
       nSolvers = 1;
   } else if (Parameters::getIntParam("wkr-strat", 1) == 3) {
      nSolvers /= 3;
   }

   if (solverType == "glucose") {
      SolverFactory::createGlucoseSolvers(nSolvers, solvers);
   } else if (solverType == "lingeling") {
      SolverFactory::createLingelingSolvers(nSolvers, solvers);
   } else if (solverType == "maple") {
      SolverFactory::createMapleSolvers(nSolvers, solvers);
   } else if (solverType == "maplesat") {
      SolverFactory::createMaplesatSolvers(nSolvers, solvers);
   } else if (solverType == "combo") {
      SolverFactory::createComboSolvers(nSolvers, solvers);
   } else {
      // MiniSat is the default choice
      SolverFactory::createMiniSatSolvers(nSolvers, solvers);
   }

   // Diversifycation
   int diversification = Parameters::getIntParam("d", 0);

   switch (diversification) {
      case 1 :
         SolverFactory::sparseDiversification(solvers);
         break;

      case 2 :
         SolverFactory::binValueDiversification(solvers);
         break;

      case 3 :
         SolverFactory::randomDiversification(solvers, 2015);
         break;

      case 4 :
         SolverFactory::nativeDiversification(solvers);
         break;

      case 5 :
         SolverFactory::sparseDiversification(solvers);
         SolverFactory::nativeDiversification(solvers);
         break;

      case 6 :
         SolverFactory::sparseRandomDiversification(solvers);
         break;

      case 7 :
         SolverFactory::sparseRandomDiversification(solvers);
         SolverFactory::nativeDiversification(solvers);
         break;

      case 0 :
         break;
   }

   if (Parameters::getIntParam("wkr-strat", 1) == 3) {
      solvers.push_back(SolverFactory::createLingelingSolver());
   }

   vector<SolverInterface *> from;
   // Start sharing threads
   switch(Parameters::getIntParam("shr-strat", 0)) {
      case 1 :
         nSharers   = 1;
         sharers    = new Sharer*[nSharers];
         sharers[0] = new Sharer(0, new SimpleSharing(), solvers, solvers);
         break;
      case 2 :
         nSharers = cpus;
         sharers  = new Sharer*[nSharers];

         for (size_t i = 0; i < nSharers; i++) {
            from.clear();
            from.push_back(solvers[i]);
            sharers[i] = new Sharer(i, new HordeSatSharing(), from,
                                    solvers);
         }
         break;
      case 3 :
         nSharers = 1;
         sharers  = new Sharer*[nSharers];
         sharers[0] = new Sharer(0, new HordeSatSharing(), solvers, solvers);
         break;

      case 4:
         nSharers = 2;
         sharers  = new Sharer*[nSharers];

         for (size_t i=0; i < nSolvers; i++) {
            from.push_back(solvers[i]);
         }
         sharers[0] = new Sharer(0, new SimpleSharing(), from, solvers);

         from.clear();
         from.push_back(solvers[nSolvers]);
         sharers[1] = new Sharer(1, new HordeSatSharing(), from, solvers);
         break;

      case 0 :
         break;
   }

   WorkingStrategy * childPF, *childCC;
   // Working strategy creation
   switch(Parameters::getIntParam("wkr-strat", 1)) {
      case 1 :
         working = new Portfolio();
         for (size_t i = 0; i < cpus; i++) {
            working->addSlave(new SequentialWorker(solvers[i]));
         }
         break;

      case 2 :
         working = new CubeAndConquer(cpus);
         working->addSlave(new SequentialWorker(solvers[0]));
         break;

      case 3 :
         working = new Portfolio();

         childPF = new Portfolio();
         for (size_t i = 0; i < nSolvers; i++) {
            childPF->addSlave(new SequentialWorker(solvers[i]));
         }
         working->addSlave(childPF);

         childCC = new CubeAndConquer(cpus - nSolvers);
         childCC->addSlave(new SequentialWorker(solvers[nSolvers]));
         working->addSlave(childCC);
         break;

      case 4 :
         working = new DivideAndConquer();
         if(Parameters::getIntParam("copy-mode",1) == 2) {
            working->addSlave(new SequentialWorker(solvers[0]));
            for(size_t i = 1; i < cpus; i++) {
	            working->addSlave(new SequentialWorker(NULL));
            }
         } else {
            for(size_t i = 0; i < cpus; i++) {
	            working->addSlave(new SequentialWorker(solvers[i]));
            }
         }
         break;

      case 0 :
         break;
   }

   // Init the management of clauses
   ClauseManager::initClauseManager();

   // Launch working
   vector<int> cube;
   working->solve(cube);

   // Wait until end or timeout
   int timeout   = Parameters::getIntParam("t", -1);
   int maxMemory = Parameters::getIntParam("max-memory", -1) * 1024 * 1024;

   while(globalEnding == false) {
      sleep(1);

      if (maxMemory > 0 && getMemoryUsed() > maxMemory) {
         cout << "c Memory used is going too large!!!!" << endl;
      }

      if (timeout > 0 && getRelativeTime() >= timeout) {
         globalEnding = true;
         working->setInterrupt();
      }
   }

   // Delete sharers
   for (int i=0; i < nSharers; i++) {
      delete sharers[i];
   }
   delete sharers;

   // Delete working strategy
   delete working;

   // Delete shared clauses
   ClauseManager::joinClauseManager();

   // Print solver stats
   //SolverFactory::printStats(solvers);

   // Print the result and the model if SAT
   log(0, "Finished\n");
   if (finalResult == SAT) {
      printf("s SATISFIABLE\n");

      if (Parameters::isSet("no-model") == false) {
         printModel(finalModel);
      }
   } else if (finalResult == UNSAT) {
      printf("s UNSATISFIABLE\n");
   } else {
      printf("s UNKNOWN\n");
   }

   return finalResult;
}
