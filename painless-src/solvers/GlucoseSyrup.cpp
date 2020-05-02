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

// Glucose includes
#include "utils/System.h"
#include "core/Dimacs.h"
#include "parallel/ParallelSolver.h"

#include "../utils/Parameters.h"
#include "../clauses/ClauseManager.h"
#include "../solvers/GlucoseSyrup.h"

using namespace Glucose;

// Macros to converte glucose literals
#define GLUE_LIT(lit) lit > 0 ? mkLit(lit - 1, false) : mkLit((-lit) - 1, true)

#define INT_LIT(lit) sign(lit) ? -(var(lit) + 1) : (var(lit) + 1)


void makeGlueVec(ClauseExchange * cls, vec<Lit> & gcls)
{
   for (size_t i = 0; i < cls->size; i++) {
      gcls.push(GLUE_LIT(cls->lits[i]));
   }
}


void glucoseExportUnary(void * issuer, Lit & l)
{
   GlucoseSyrup * gs = (GlucoseSyrup*)issuer;

   ClauseExchange * ncls = ClauseManager::allocClause(1);

   ncls->from    = gs->id;
   ncls->lbd     = 1;
   ncls->lits[0] = INT_LIT(l);

   gs->clausesToExport.addClause(ncls); 
}

void glucoseExportClause(void * issuer, Clause & cls)
{
   GlucoseSyrup * gs = (GlucoseSyrup *)issuer;

   if (cls.lbd() > gs->glueLimit || gs->glueLimit <= 0)
      return;

   ClauseExchange * ncls = ClauseManager::allocClause(cls.size());

   ncls->from = gs->id;
   ncls->lbd  = cls.lbd();

   for (int i = 0; i < cls.size(); i++) {
      ncls->lits[i] = INT_LIT(cls[i]);
   }

   gs->clausesToExport.addClause(ncls); 
}

Lit glucoseImportUnary(void * issuer)
{
   GlucoseSyrup * gs = (GlucoseSyrup*)issuer;

   Lit l = lit_Undef;

   ClauseExchange * cls = NULL;

   if (gs->unitsToImport.getClause(&cls) == false)
      return l;

   l = GLUE_LIT(cls->lits[0]);

   ClauseManager::releaseClause(cls);

   return l;
}

bool glucoseImportClause(void * issuer, int * from, vec<Lit> & gcls)
{
   GlucoseSyrup * gs = (GlucoseSyrup*)issuer;

   ClauseExchange * cls = NULL;

   if (gs->clausesToImport.getClause(&cls) == false)
      return false;

   makeGlueVec(cls, gcls);

   *from = cls->from;

   ClauseManager::releaseClause(cls);

   return true;
}

GlucoseSyrup::GlucoseSyrup(int id) : SolverInterface(id, GLUCOSE)
{
   glueLimit = Parameters::getIntParam("lbd-limit", 100);

   solver = new ParallelSolver(id);

   switch(Parameters::getIntParam("split-heur",1)) {
      case 2:
         solver->useFlip=true;
         break;
      case 3:
         solver->usePR=true;
         break;
      default:;
   }

   solver->exportUnary  = glucoseExportUnary;
   solver->exportClause = glucoseExportClause;
   solver->importUnary  = glucoseImportUnary;
   solver->importClause = glucoseImportClause;
   solver->issuer       = this;
}

GlucoseSyrup::GlucoseSyrup(const GlucoseSyrup & other, int id) :
   SolverInterface(id, GLUCOSE)
{
   glueLimit = Parameters::getIntParam("lbd-limit", 100);

   solver = new ParallelSolver(*(other.solver), id);

   switch(Parameters::getIntParam("split-heur",1)) {
      case 2:
         solver->useFlip=true;
         break;
      case 3:
         solver->usePR=true;
         break;
      default:;
   }
   
   solver->exportUnary  = glucoseExportUnary;
   solver->exportClause = glucoseExportClause;
   solver->importUnary  = glucoseImportUnary;
   solver->importClause = glucoseImportClause;
   solver->issuer       = this;
}

GlucoseSyrup::~GlucoseSyrup()
{
  vector<ClauseExchange *> tmp;
  
  unitsToImport.getClauses(tmp);
  clausesToImport.getClauses(tmp);
  clausesToExport.getClauses(tmp);
  clausesToAdd.getClauses(tmp);
  
  for(size_t i=0;i<tmp.size();i++){
    if(tmp[i]==NULL){ printf("NAON :'(\n");abort();}
    ClauseManager::releaseClause(tmp[i]);
    }
   delete solver;
}

bool
GlucoseSyrup::loadFormula(const char * filename)
{
   gzFile in = gzopen(filename, "rb");

   parse_DIMACS(in, *solver);

   gzclose(in);

   return true;
}

//Get the number of variables of the formula
int
GlucoseSyrup::getVariablesCount()
{
   return solver->nVars();
}

// Get a variable suitable for search splitting
vector<int>
GlucoseSyrup::getDivisionVariables(int k) // NOTE: ignoring values of k > 1 for now
{
  Lit res;

  switch(Parameters::getIntParam("split-heur",1)) {
  case 2:
    res=solver->pickBranchLitUsingFlipActivity();
    break;
  case 3:
    res=solver->pickBranchLitUsingPropagationRate();
    break;
  default:
    res=solver->pickBranchLit();
  }

  return vector<int>(1, INT_LIT(res));
}

// Set initial phase for a given variable
void
GlucoseSyrup::setPhase(const int var, const bool phase)
{
   solver->setPolarity(var-1, phase);
}

//Bump activity for a given variable
void
GlucoseSyrup::bumpVariableActivity(const int var, const int times)
{
   for(int i = 0; i < times; i++) {
      solver->varBumpActivity(var-1);
   } 
}

// Interrupt the SAT solving, so it can be started again with new assumptions
void
GlucoseSyrup::setSolverInterrupt()
{
   solver->interrupt();
}

// Diversify the solver
void
GlucoseSyrup::diversify(int id)
{
   int idMod = id ? id <= 9 : id % 8;

   switch(idMod) {
      case 1 :
         solver->var_decay     = 0.94;
         solver->max_var_decay = 0.96;
         solver->firstReduceDB = 600;
         break;

      case 2 :
         solver->var_decay     = 0.90;
         solver->max_var_decay = 0.97;
         solver->firstReduceDB = 500;
         break;

      case 3 :
         solver->var_decay     = 0.85;
         solver->max_var_decay = 0.93;
         solver->firstReduceDB = 400;
         break;

      case 4 :
         solver->var_decay       = 0.95;
         solver->max_var_decay   = 0.95;
         solver->firstReduceDB   = 4000;
         solver->sizeLBDQueue    = 100;
         solver->K               = 0.7;
         solver->incReduceDB     = 500;
         solver->lbdQueue.growTo (100);
         break;

      case 5 :
         solver->var_decay     = 0.93;
         solver->max_var_decay = 0.96;
         solver->firstReduceDB = 100;
         solver->incReduceDB   = 500;
         break;

      case 6 :
         solver->var_decay     = 0.75;
         solver->max_var_decay = 0.94;
         solver->firstReduceDB = 2000;
         break;

      case 7 :
         solver->var_decay     = 0.94;
         solver->max_var_decay = 0.96;
         solver->firstReduceDB = 800;
         break;

      case 8 :
         solver->reduceOnSize = true;
         break;

      case 9 :
         solver->reduceOnSize     = true;
         solver->reduceOnSizeSize = 14;
         break;

      case 0 :
      default:
         break;
   }

   if (id > 9) {
      int noiseFactor       = id / 8;
      double noisevar_decay = 0.005 + noiseFactor * 0.006;
      int noiseReduceDB     = 50 + noiseFactor * 25;

      solver->var_decay     += noisevar_decay;
      solver->firstReduceDB +=noiseReduceDB;
   }
}

void
GlucoseSyrup::unsetSolverInterrupt()
{
   solver->clearInterrupt();
}

// Solve the formula with a given set of assumptions
// return 10 for SAT, 20 for UNSAT, 0 for UNKNOWN
SatResult
GlucoseSyrup::solve(const vector<int> & cube, int nConflicts)
{
   unsetSolverInterrupt();

   vector<ClauseExchange *> tmp;
   clausesToAdd.getClauses(tmp);

   for (size_t i = 0; i < tmp.size(); i++) {
      vec<Lit> gcls;
      makeGlueVec(tmp[i], gcls);

      ClauseManager::releaseClause(tmp[i]);

      if (solver->addClause(gcls) == false) {
         printf("unsat when adding cls\n");
         return UNSAT;
      }
   }

   vec<Lit> gAssumptions;
   for (int i = 0; i < cube.size(); i++) {
     Lit l = GLUE_LIT(cube[i]);
     if (!solver->isEliminated(var(l))) {
         gAssumptions.push(l);
      }
   }

   if ( nConflicts != -1 )
       solver->setConfBudget(nConflicts);

   lbool res = solver->solveLimited(gAssumptions);

   if (res == l_True)
      return SAT;

   if (res == l_False)
      return UNSAT;

   return UNKNOWN;
}

void
GlucoseSyrup::addClause(ClauseExchange * clause)
{
   clausesToAdd.addClause(clause);

   setSolverInterrupt();
}

void
GlucoseSyrup::addLearnedClause(ClauseExchange * clause)
{
   if (clause->size == 1) {
      unitsToImport.addClause(clause);
   } else {
      clausesToImport.addClause(clause);
   }
}

void
GlucoseSyrup::addClauses(const vector<ClauseExchange *> & clauses)
{
   clausesToAdd.addClauses(clauses);

   setSolverInterrupt();
}

void
GlucoseSyrup::addInitialClauses(const vector<ClauseExchange *> & clauses)
{
   for (size_t ind = 0; ind < clauses.size(); ind++) {
      vec<Lit> mcls;

      for (size_t i = 0; i < clauses[ind]->size; i++) {
         int lit = clauses[ind]->lits[i];
         int var = abs(lit);

         while (solver->nVars() < var) {
            solver->newVar();
         }

         mcls.push(GLUE_LIT(lit));
      }

      if (!solver->addClause(mcls)) {
         printf("unsat when adding initial cls\n");
      }
   }
}

void
GlucoseSyrup::addLearnedClauses(const vector<ClauseExchange *> & clauses)
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
GlucoseSyrup::getLearnedClauses(vector<ClauseExchange *> & clauses)
{
   clausesToExport.getClauses(clauses);
}

void
GlucoseSyrup::increaseClauseProduction()
{
   glueLimit++;
}

void
GlucoseSyrup::decreaseClauseProduction()
{
   glueLimit--;
}

SolvingStatistics
GlucoseSyrup::getStatistics()
{
   SolvingStatistics stats;

   stats.conflicts    = solver->conflicts;
   stats.propagations = solver->propagations;
   stats.restarts     = solver->starts;
   stats.decisions    = solver->decisions;
   stats.memPeak      = memUsedPeak();

   return stats;
}

std::vector<int>
GlucoseSyrup::getModel()
{
   std::vector<int> model;

   for (int i = 0; i < solver->nVars(); i++) {
      if (solver->model[i] != l_Undef) {
         int lit = solver->model[i] == l_True ? i + 1 : -(i + 1);

         model.push_back(lit);
      }
   }

   return model;
}

ClauseExchange *
GlucoseSyrup::getFinalAnalysis()
{
   ClauseExchange * out_cls =
      ClauseManager::allocClause(solver->conflict.size());

   out_cls->lbd = 2;

   for (int i = 0; i < solver->conflict.size(); i++) {
      out_cls->lits[i] = INT_LIT(solver->conflict[i]);
   }

   return out_cls;
}
