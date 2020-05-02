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

#include "../clauses/ClauseManager.h"
#include "../utils/Parameters.h"
#include "../solvers/Lingeling.h"

// Include lingeling
extern "C" {
	#include "lglib.h"
}

#include <ctype.h>

int termCallback(void * solverPtr)
{
	Lingeling * lp = (Lingeling *)solverPtr;
	return lp->stopSolver;
}

void produceUnit(void * sp, int lit)
{
   Lingeling * lp = (Lingeling *)sp;
   
   // Create new clause
   ClauseExchange * ncls = ClauseManager::allocClause(1);

   ncls->lits[0] = lit;
   ncls->from    = lp->id;

   // Add it to the buffer for export
   lp->clausesToExport.addClause(ncls);
}

void produce(void * sp, int * cls, int glue)
{
   // If unit clause, call produceUnit
   if (cls[1] == 0) {
      produceUnit(sp, cls[0]);
      return;
   }

   Lingeling * lp = (Lingeling *)sp;

   // If the clause has not a good lbd
   if (glue > lp->glueLimit || lp->glueLimit < 0)
      return;

   // Create new clause
   int size = 0;
   while (cls[size] != 0) {
      size++;
   }

   ClauseExchange * ncls = ClauseManager::allocClause(size);

   memcpy(ncls->lits, cls, sizeof(int)*size);

   ncls->lbd  = glue;
   ncls->from = lp->id;

   // Add it to the buffer for export
   lp->clausesToExport.addClause(ncls);
}

void consumeUnits(void * sp, int ** start, int ** end)
{
   Lingeling* lp = (Lingeling*)sp;
   
   vector<ClauseExchange *> tmp;

   lp->unitsToImport.getClauses(tmp);

	if (tmp.empty()) {
		*start = lp->unitsBuffer;
		*end   = lp->unitsBuffer;
		return;
	}

	if (tmp.size() >= lp->unitsBufferSize) {
		lp->unitsBufferSize = 1.6 * tmp.size();
		lp->unitsBuffer     = (int *)realloc((void *)lp->unitsBuffer,
                                           lp->unitsBufferSize * sizeof(int));
	}

	for (size_t i = 0; i < tmp.size(); i++) {
		lp->unitsBuffer[i] = tmp[i]->lits[0];
      ClauseManager::releaseClause(tmp[i]);
	}

	*start = lp->unitsBuffer;
	*end   = *start + tmp.size();
}

void consumeCls(void* sp, int ** clause, int * glue)
{
   Lingeling* lp = (Lingeling*)sp;

   ClauseExchange * cls = NULL;
   
   if (lp->clausesToImport.getClause(&cls) == false) {
		*clause = NULL;
		return;
	}

   if (cls->size+1 >= lp->clsBufferSize) {
      lp->clsBufferSize = 1.6 * cls->size;
      lp->clsBuffer     = (int *)realloc((void *)lp->clsBuffer,
                                         lp->clsBufferSize * sizeof(int));
   }
	
   *glue = cls->lbd;

   memcpy(lp->clsBuffer, cls->lits, sizeof(int)*cls->size);
   lp->clsBuffer[cls->size] = 0;

   ClauseManager::releaseClause(cls);

   *clause = lp->clsBuffer;
}

Lingeling::Lingeling(int id) : SolverInterface(id, LINGELING)
{
   solver = lglinit();

   // BCA has to be disabled for valid clause sharing (or freeze all literals)
   lglsetopt(solver, "druplig", 0);
   lglsetopt(solver, "bca", 0);
   lglsetopt(solver, "profile", 0);

   stopSolver = 0;

   glueLimit = Parameters::getIntParam("lbd-limit", 2);
   
   unitsBufferSize = clsBufferSize = 100;
   unitsBuffer     = (int*) malloc(unitsBufferSize * sizeof(int));
   clsBuffer       = (int*) malloc(clsBufferSize * sizeof(int));
	
   lglsetproducecls  (solver, produce, this);
   lglsetproduceunit (solver, produceUnit, this);
   lglsetconsumeunits(solver, consumeUnits, this);
   lglsetconsumecls  (solver, consumeCls, this);
   lglseterm         (solver, termCallback, this);
}

Lingeling::Lingeling(const Lingeling & other, int id) :
   SolverInterface(id, LINGELING)
{
   solver = lglclone(other.solver);

   // BCA has to be disabled for valid clause sharing (or freeze all literals)
   lglsetopt(solver, "druplig", 0);
   lglsetopt(solver, "bca", 0);
   lglsetopt(solver, "profile", 0);

   stopSolver = 0;

   glueLimit = Parameters::getIntParam("lbd-limit", 2);
   
   unitsBufferSize = clsBufferSize = 100;
   unitsBuffer     = (int*) malloc(unitsBufferSize * sizeof(int));
   clsBuffer       = (int*) malloc(clsBufferSize * sizeof(int));
	
   lglsetproducecls  (solver, produce, this);
   lglsetproduceunit (solver, produceUnit, this);
   lglsetconsumeunits(solver, consumeUnits, this);
   lglsetconsumecls  (solver, consumeCls, this);
   lglseterm         (solver, termCallback, this);
}

Lingeling::~Lingeling()
{
   lglrelease(solver);

   free(unitsBuffer);
   free(clsBuffer);
}

bool
Lingeling::loadFormula(const char* filename)
{
   vector<SolverInterface*> solversToManage;

   solversToManage.push_back(this);

   bool ret = loadFormulaToSolvers(solversToManage, filename);

   lglsimp(solver, 10);

   return ret;
}

int
Lingeling::getVariablesCount()
{
   return lglnvars(solver);
}

// Get a variable suitable for search splitting
vector<int>
Lingeling::getDivisionVariables(int k)
{
   lglsimp(solver, 1);

   int oldjwhred = lglgetopt (solver, "jwhred");

   int lit = lglookahead(solver);

   lglsetopt(solver, "jwhred", oldjwhred);

   return vector<int>(1, lit);
}

// Set initial phase for a given variable
void
Lingeling::setPhase(const int var, const bool phase)
{
	lglsetphase(solver, phase ? var : -var);
}

// Bump activity for a given variable
void
Lingeling::bumpVariableActivity(const int var, const int times)
{
   for (int i = times; i < times; i++) {
      // TODO: is it good ?
      lglbumpdlit(solver, var);
   }
}

// Interrupt the SAT solving, so it can be started again with new assumptions
void
Lingeling::setSolverInterrupt()
{
   stopSolver = 1;
}

void
Lingeling::unsetSolverInterrupt()
{
   stopSolver = 0;
}

// Solve the formula with a given set of assumptions
// return 10 for SAT, 20 for UNSAT, 0 for UNKNOWN
SatResult
Lingeling::solve(const vector<int> & cube, int nConflicts)
{
   // add the clauses
   vector<ClauseExchange *> tmp;
   clausesToAdd.getClauses(tmp);

   for (size_t i = 0; i < tmp.size(); i++) {
      for (size_t j = 0; j < tmp[i]->size; j++) {
         lgladd(solver, tmp[i]->lits[j]);
      }

      lgladd(solver, 0);

      ClauseManager::releaseClause(tmp[i]);
   }

   // set the assumptions
   for (size_t i = 0; i < cube.size(); i++) {
      // freezing problems ???
      if (lglusable(solver, cube[i])) {
         lglassume(solver, cube[i]);
      }
   }

   if (nConflicts) {
      cerr << "n_conflicts limitation solving is not yet implemented for" \
               "Lingeling solver" << endl;
      exit(1);
   }

   // Simplify the problem
   int res = lglsimp(solver, 0);

   switch (res) {
      case LGL_SATISFIABLE:
         return SAT;
      case LGL_UNSATISFIABLE:
         return UNSAT;
   }

   // Solve the problem
   res = lglsat(solver);

   switch (res) {
      case LGL_SATISFIABLE:
         return SAT;
      case LGL_UNSATISFIABLE:
         return UNSAT;
   }

   return UNKNOWN;
}

// Add a permanent clause to the formula
void
Lingeling::addClause(ClauseExchange * clause)
{
   clausesToAdd.addClause(clause);
}

void
Lingeling::addClauses(const vector<ClauseExchange *> & clauses)
{
   clausesToAdd.addClauses(clauses);
}

void
Lingeling::addInitialClauses(const vector<ClauseExchange *> & clauses)
{
   for (size_t i = 0; i < clauses.size(); i++) {
      for (size_t j = 0; j < clauses[i]->size; j++) {
         lgladd(solver, clauses[i]->lits[j]);
      }

      lgladd(solver, 0);
   }
}

// Add a learned clause to the formula
void
Lingeling::addLearnedClause(ClauseExchange * clause)
{
   if (clause->size == 1) {
      unitsToImport.addClause(clause);
   } else {
      clausesToImport.addClause(clause);
   }
}

void
Lingeling::addLearnedClauses(const vector<ClauseExchange *> & clauses)
{
   for (size_t i = 0; i < clauses.size(); i++) {
      if (clauses[i]->size == 1) {
         unitsToImport.addClause(clauses[i]);
      } else {
         clausesToImport.addClause(clauses[i]);
      }
   }
}

void
Lingeling::increaseClauseProduction()
{
   glueLimit++;
}

void
Lingeling::decreaseClauseProduction()
{
   glueLimit--;
}

void
Lingeling::getLearnedClauses(vector<ClauseExchange *> & clauses)
{
   return clausesToExport.getClauses(clauses);
}

SolvingStatistics
Lingeling::getStatistics()
{
   SolvingStatistics stats;

   stats.conflicts    = lglgetconfs(solver);
   stats.decisions    = lglgetdecs(solver);
   stats.propagations = lglgetprops(solver);
   stats.memPeak      = lglmaxmb(solver);

   return stats;
}

void
Lingeling::diversify(int id)
{
   // This method is copied from Plingeling
   lglsetopt(solver, "seed", id);
   lglsetopt (solver, "classify", 0);

   switch (id % 13) {
      case 0 :
      default :
         break;

      case 1 :
         lglsetopt(solver, "plain", 1),
         lglsetopt(solver, "decompose", 1);
         break;

      case 2 :
         lglsetopt(solver, "restartint", 1000);
         break;

      case 3 :
         lglsetopt(solver, "elmresched", 7);
         break;

      case 4 :
         lglsetopt(solver, "scincincmin", 250);
         break;

      case 5 :
         lglsetopt(solver, "block", 0);
         lglsetopt(solver, "cce", 0);
         break;

      case 6 :
         lglsetopt(solver, "scincinc", 50);
         break;

      case 7 :
         lglsetopt(solver, "phase", -1);
         break;

      case 8 :
         lglsetopt(solver, "phase", 1);
         break;

      case 9 :
         lglsetopt(solver, "sweeprtc", 1);
         break;

      case 10 :
         lglsetopt(solver, "restartint", 100);
         break;

      case 11 :
         lglsetopt(solver, "reduceinit", 10000);
         lglsetopt(solver, "reducefixed", 1);
         break;

      case 12 :
         lglsetopt(solver, "restartint", 4);
         break;
   }
   //switch (id % 16) {
   //   case 1:
   //      lglsetopt (solver, "plain", 1);
   //      break;

   //   case 2 :
   //      lglsetopt(solver, "agilelim", 100);
   //      break;

   //   case 3 :
   //      lglsetopt(solver, "block", 0);
   //      lglsetopt(solver, "cce", 0);
   //      break;

   //   case 4 :
   //      lglsetopt(solver, "bias", -1);
   //      break;

   //   case 5 :
   //      lglsetopt(solver, "acts", 0);
   //      break;

   //   case 6 :
   //      lglsetopt(solver, "phase", 1);
   //      break;

   //   case 7 :
   //      lglsetopt(solver, "acts", 1);
   //      break;

   //   case 8 :
   //      lglsetopt(solver, "bias", 1);
   //      break;

   //   case 9 :
   //      lglsetopt(solver, "wait", 0);
   //      lglsetopt(solver, "blkrtc", 1);
   //      lglsetopt(solver, "elmrtc", 1);
   //      break;

   //   case 10 :
   //      lglsetopt(solver, "phase", -1);
   //      break;

   //   case 11 :
   //      lglsetopt(solver, "prbsimplertc", 1);
   //      break;

   //   case 12 :
   //      lglsetopt(solver, "gluescale", 1);
   //      break;

   //   case 13 :
   //      lglsetopt(solver, "gluescale", 3);
   //      break;

   //   case 14 :
   //      lglsetopt(solver, "move", 1);
   //      break;

   //   case 15 :
   //      lglsetopt(solver, "flipping", 1);
   //      break;

   //   case 0 :
   //   default :
   //      break;
   //}
}

vector<int>
Lingeling::getModel()
{
   vector<int> model;
   for (int i = 1; i <= lglmaxvar(solver); i++) {
	   int lit = (lglderef(solver, i) > 0) ? i : -i;
	   model.push_back(lit);
   }

   return model;
}

ClauseExchange *
Lingeling::getFinalAnalysis()
{
   cerr << "The function getFinalAnalysis is not implemented for the " \
            "Lingeling solver" << endl;
   exit(1);

   return NULL;
}
