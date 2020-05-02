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

// Maple includes
#include "mapleCOMSPS/utils/System.h"
#include "mapleCOMSPS/core/Dimacs.h"
#include "mapleCOMSPS/simp/SimpSolver.h"

#include "../community/GraphAlgorithms.h"
#include "../utils/Logger.h"
#include "../utils/System.h"
#include "../utils/Parameters.h"
#include "../clauses/ClauseManager.h"
#include "../solvers/Maple.h"

using namespace MapleCOMSPS;

// Macros for minisat literal representation conversion
#define MINI_LIT(lit) lit > 0 ? mkLit(lit - 1, false) : mkLit((-lit) - 1, true)

#define INT_LIT(lit) sign(lit) ? -(var(lit) + 1) : (var(lit) + 1)


static void makeMiniVec(ClauseExchange * cls, vec<Lit> & mcls)
{
   for (size_t i = 0; i < cls->size; i++) {
      mcls.push(MINI_LIT(cls->lits[i]));
   }
}


void mapleExportClause(void * issuer, int lbd, vec<Lit> & cls)
{
	Maple * mp = (Maple*)issuer;

	if (lbd > mp->lbdLimit)
		return;

	ClauseExchange * ncls = ClauseManager::allocClause(cls.size());

   ncls->lbd = lbd;

	for (int i = 0; i < cls.size(); i++) {
		ncls->lits[i] = INT_LIT(cls[i]);
	}

   ncls->from = mp->id;

   mp->clausesToExport.addClause(ncls);
}

Lit mapleImportUnit(void * issuer)
{
   Maple * mp = (Maple*)issuer;

   Lit l = lit_Undef;

   ClauseExchange * cls = NULL;

   if (mp->unitsToImport.getClause(&cls) == false)
      return l;

   l = MINI_LIT(cls->lits[0]);

   ClauseManager::releaseClause(cls);

   return l;
}

bool mapleImportClause(void * issuer, int * lbd, vec<Lit> & mcls)
{
   Maple* mp = (Maple*)issuer;

   ClauseExchange * cls = NULL;

   if (mp->clausesToImport.getClause(&cls) == false)
      return false;

   makeMiniVec(cls, mcls);

   *lbd = cls->lbd;

   ClauseManager::releaseClause(cls);

   return true;
}

Maple::Maple(int id) : SolverInterface(id, MAPLE)
{
   lbdLimit = Parameters::getIntParam("lbd-limit", 2);

   solver = new SimpSolver();

   switch(Parameters::getIntParam("split-heur",1)) {
	   case 2:
         solver->useFlip=true;
         break;
      case 3:
         solver->usePR=true;
         break;
      case 9:
         solver->useFlip = true;
         solver->usePR = true;
         break;
      case 10:
         solver->useFlip = true;
         solver->usePR = true;
         break;
      default:;
   }

   solver->exportClauseCallback = mapleExportClause;
   solver->importUnitCallback   = mapleImportUnit;
   solver->importClauseCallback = mapleImportClause;
   solver->issuer               = this;
}

Maple::Maple(const Maple & other, int id) : SolverInterface(id, MAPLE)
{
   lbdLimit = Parameters::getIntParam("lbd-limit", 2);

   solver = new SimpSolver(*(other.solver));

   switch(Parameters::getIntParam("split-heur",1)) {
      case 2:
         solver->useFlip=true;
         break;
      case 3:
         solver->usePR=true;
         break;
      case 9:
         solver->useFlip = true;
         solver->usePR = true;
         break;
      case 10:
         solver->useFlip = true;
         solver->usePR = true;
         break;
      default:;
   }
  
   solver->exportClauseCallback = mapleExportClause;
   solver->importUnitCallback   = mapleImportUnit;
   solver->importClauseCallback = mapleImportClause;
   solver->issuer               = this;
}

Maple::~Maple()
{
	delete solver;
}

bool
Maple::loadFormula(const char * filename)
{
    gzFile in = gzopen(filename, "rb");

    parse_DIMACS(in, *solver);

    gzclose(in);

    bool use_gaussian = (Parameters::isSet("gaussian") == true);
    solver->eliminate(false, use_gaussian);

    if (id % 3 == 0)
    {
        solver->init_bayesian();
        solver->bayesian();
    }
    else if (id % 3 == 1)
    {
        solver->jeroslow_wang_init(true, true);
    }
//    if ( Parameters::getIntParam("split-heur", 1) == 10 )
//    {
//        solver->init_bayesian();
//        solver->bayesian();
//    }
//
//    bool jw_act = (Parameters::getIntParam("init-act", 0) == 1);
//    bool jw_pol = (Parameters::getIntParam("init-pol", 0) == 1);
//    if ( jw_act || jw_pol )
//        solver->jeroslow_wang_init(jw_act, jw_pol);
    
    

    return true;
}

//Get the number of variables of the formula
int
Maple::getVariablesCount()
{
	return solver->nVars();
}

// Get a variable suitable for search splitting
vector<int>
Maple::getDivisionVariables(int k) // NOTE: ignoring values of k > 1 for now
{
   Lit res;
   vector<map<int, double> > al;
   Graph vig;
   double * page_rank, * ordered_page_rank;
   vector<float> kcore, ordered_kcore;
   vector<int> not_used;
   int nVars, next = 1;

   switch(Parameters::getIntParam("split-heur",1)) {
      case 2:
        res=solver->pickBranchLitUsingFlipActivity();
        break;

      case 3:
        res=solver->pickBranchLitUsingPropagationRate();
        break;

      case 4:
        return vector<int>(1, (rand() % getVariablesCount()) + 1);

      case 5:
         al = solver->createVIG();
         vig = Graph(al);
         nVars = al.size();

         page_rank = new double[nVars];
         pageRank(vig, page_rank);

         ordered_page_rank = new double[nVars];
         copy(page_rank, page_rank + nVars, ordered_page_rank);
         sort(ordered_page_rank, ordered_page_rank + nVars);

         for (int i = nVars - 1; i >= 0; i--) {
           next = distance(page_rank, find(page_rank, page_rank + nVars,
                  ordered_page_rank[i]));
            if (solver->decision[next]) {
                next++;
                break;
            }
         }
         return vector<int>(1, next);

      case 6:
         al = solver->createVIG();
         vig = Graph(al);
         nVars = al.size();

         page_rank = new double[nVars];
         weightedPageRank(vig, page_rank);

         ordered_page_rank = new double[nVars];
         copy(page_rank, page_rank + nVars, ordered_page_rank);
         sort(ordered_page_rank, ordered_page_rank + nVars);

         for (int i = nVars - 1; i >= 0; i--) {
            next = distance(page_rank, find(page_rank, page_rank + nVars,
                   ordered_page_rank[i]));
            if (solver->decision[next]) {
                next++;
                break;
            }
         }
         return vector<int>(1, next);

      case 7:
         al = solver->createVIG();
         vig = Graph(al);
         nVars = al.size();
         
         kcore = weightedKCore(vig, nVars, not_used);

         ordered_kcore.resize(nVars, 0.);
         copy(kcore.begin(), kcore.end(), ordered_kcore.begin());
         sort(ordered_kcore.begin(), ordered_kcore.end());

         for (int i = nVars - 1; i >= 0; i--) {
            next = distance(kcore.begin(), find(kcore.begin(), kcore.end(), ordered_kcore[i]));
            if (solver->decision[next]) {
                next++;
                break;
            }
          }
          return vector<int>(1, next);

      case 9:
          return solver->pickSplittingVariables();

      case 10: // bayesian splitting
          res = solver->pickBranchLitUsingBayesian();
          break;
   
      case 11: // LRB
         return solver->pickSplittingVariablesUsingLRBLit(1);

      default:
        res = solver->pickBranchLit();
   }

   return vector<int>(1, INT_LIT(res));
}

// Set initial phase for a given variable
void
Maple::setPhase(const int var, const bool phase)
{
   solver->setPolarity(var - 1, phase ? true : false);
}

// Bump activity for a given variable
void
Maple::bumpVariableActivity(const int var, const int times)
{
   for(int i = 0; i < times; i++) {
      solver->varBumpActivity(var - 1, 1);
      //TODO: work only for VSIDS
   }
}

// Interrupt the SAT solving, so it can be started again with new assumptions
void
Maple::setSolverInterrupt()
{
   stopSolver = true;

	solver->interrupt();
}

void
Maple::unsetSolverInterrupt()
{
   stopSolver = false;

	solver->clearInterrupt();
}

// Diversify the solver
void
Maple::diversify(int id)
{
	//solver->random_seed = (double)id;
   if (id % 2) {
      solver->VSIDS = true;
   } else {
      solver->VSIDS = false;
   }
}

// Solve the formula with a given set of assumptions
// return 10 for SAT, 20 for UNSAT, 0 for UNKNOWN
SatResult
Maple::solve(const vector<int> & cube, int nConflicts)
{
   unsetSolverInterrupt();

   vector<ClauseExchange *> tmp;
   
   tmp.clear();
   clausesToAdd.getClauses(tmp);

   for (size_t ind = 0; ind < tmp.size(); ind++) {
      vec<Lit> mcls;
      makeMiniVec(tmp[ind], mcls);

      ClauseManager::releaseClause(tmp[ind]);

      if (solver->addClause(mcls) == false) {
         printf("c unsat when adding cls\n");
         return UNSAT;
      }
   }

   vec<Lit> miniAssumptions;
   for (size_t ind = 0; ind < cube.size(); ind++) {
     Lit l= MINI_LIT(cube[ind]);
     if(!solver->isEliminated(var(l))){
       miniAssumptions.push(l);
     }
   }

   if ( nConflicts != -1 )
       solver->setConfBudget(nConflicts);

   lbool res = solver->solveLimited(miniAssumptions);

   if (res == l_True)
      return SAT;

   if (res == l_False)
      return UNSAT;

   return UNKNOWN;
}

void
Maple::addClause(ClauseExchange * clause)
{
   clausesToAdd.addClause(clause);

   setSolverInterrupt();
}

void
Maple::addLearnedClause(ClauseExchange * clause)
{
   if (clause->size == 1) {
      unitsToImport.addClause(clause);
   } else {
      clausesToImport.addClause(clause);
   }
}

void
Maple::addClauses(const vector<ClauseExchange *> & clauses)
{
   clausesToAdd.addClauses(clauses);

   setSolverInterrupt();
}

void
Maple::addInitialClauses(const vector<ClauseExchange *> & clauses)
{
   for (size_t ind = 0; ind < clauses.size(); ind++) {
      vec<Lit> mcls;

      for (size_t i = 0; i < clauses[ind]->size; i++) {
         int lit = clauses[ind]->lits[i];
         int var = abs(lit);

         while (solver->nVars() < var) {
            solver->newVar();
         }

         mcls.push(MINI_LIT(lit));
      }

      if (solver->addClause(mcls) == false) {
         printf("c unsat when adding initial cls\n");
      }
   }
}

void
Maple::addLearnedClauses(const vector<ClauseExchange *> & clauses)
{
   for (size_t i = 0; i < clauses.size(); i++) {
      addLearnedClause(clauses[i]);
   }
}

void
Maple::getLearnedClauses(vector<ClauseExchange *> & clauses)
{
   clausesToExport.getClauses(clauses);
}

void
Maple::increaseClauseProduction()
{
   lbdLimit++;
}

void
Maple::decreaseClauseProduction()
{
   if (lbdLimit > 2) {
      lbdLimit--;
   }
}

SolvingStatistics
Maple::getStatistics()
{
   SolvingStatistics stats;

   stats.conflicts    = solver->conflicts;
   stats.propagations = solver->propagations;
   stats.restarts     = solver->starts;
   stats.decisions    = solver->decisions;
   stats.memPeak      = memUsedPeak();

   return stats;
}

vector<int>
Maple::getModel()
{
   vector<int> model;

   for (int i = 0; i < solver->nVars(); i++) {
      if (solver->model[i] != l_Undef) {
         int lit = solver->model[i] == l_True ? i + 1 : -(i + 1);
         model.push_back(lit);
      }
   }

   return model;
}

ClauseExchange *
Maple::getFinalAnalysis()
{
   ClauseExchange * out_cls =
      ClauseManager::allocClause(solver->conflict.size());

   out_cls->lbd = 2;

   for (int i = 0; i < solver->conflict.size(); i++) {
      out_cls->lits[i] = INT_LIT(solver->conflict[i]);
   }

   return out_cls;
}
