// Maplesat headers
#include "maplesat/utils/System.h"
#include "maplesat/core/Dimacs.h"
#include "maplesat/simp/SimpSolver.h"

//#include "../community/GraphAlgorithms.h"
//#include "../utils/Logger.h"
//#include "../utils/System.h"
#include "../utils/Parameters.h"
#include "../clauses/ClauseManager.h"
#include "../solvers/Maplesat.h"

using namespace MapleSAT;

// Macros for minisat literal representation conversion
#define MINI_LIT(lit) lit > 0 ? mkLit(lit - 1, false) : mkLit((-lit) - 1, true)

#define INT_LIT(lit) sign(lit) ? -(var(lit) + 1) : (var(lit) + 1)


static void makeMiniVec(ClauseExchange * cls, vec<Lit> & mcls)
{
    for (size_t i = 0; i < cls->size; i++) {
        mcls.push(MINI_LIT(cls->lits[i]));
    }
}


void maplesatExportClause(void * issuer, int lbd, vec<Lit> & cls)
{
    Maplesat * mp = (Maplesat*)issuer;

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

Lit maplesatImportUnit(void * issuer)
{
    Maplesat * mp = (Maplesat*)issuer;

    Lit l = lit_Undef;

    ClauseExchange * cls = NULL;

    if (mp->unitsToImport.getClause(&cls) == false)
        return l;

    l = MINI_LIT(cls->lits[0]);

    ClauseManager::releaseClause(cls);

    return l;
}

bool maplesatImportClause(void * issuer, int * lbd, vec<Lit> & mcls)
{
    Maplesat* mp = (Maplesat*)issuer;

    ClauseExchange * cls = NULL;

    if (mp->clausesToImport.getClause(&cls) == false)
        return false;

    makeMiniVec(cls, mcls);

    *lbd = cls->lbd;

    ClauseManager::releaseClause(cls);

    return true;
}

Maplesat::Maplesat(int id) : SolverInterface(id, MAPLESAT)
{
    lbdLimit = Parameters::getIntParam("lbd-limit", 2);

    solver = new SimpSolver();

    solver->exportClauseCallback = maplesatExportClause;
    solver->importUnitCallback   = maplesatImportUnit;
    solver->importClauseCallback = maplesatImportClause;
    solver->issuer               = this;
}

Maplesat::Maplesat(const Maplesat & other, int id) : SolverInterface(id, MAPLESAT)
{
   lbdLimit = Parameters::getIntParam("lbd-limit", 2);

   solver = new SimpSolver(*(other.solver));
  
   solver->exportClauseCallback = maplesatExportClause;
   solver->importUnitCallback   = maplesatImportUnit;
   solver->importClauseCallback = maplesatImportClause;
   solver->issuer               = this;
}


Maplesat::~Maplesat()
{
    delete solver;
}

bool
Maplesat::loadFormula(const char * filename)
{
    gzFile in = gzopen(filename, "rb");

    parse_DIMACS(in, *solver);

    gzclose(in);

    solver->eliminate();

    return true;
}

//Get the number of variables of the formula
int
Maplesat::getVariablesCount()
{
    return solver->nVars();
}

// Get a variable suitable for search splitting
vector<int>
Maplesat::getDivisionVariables(int k)
{
    //return solver->pickSplittingVariables(k);
    
    Lit res = solver->pickBranchLit();
    return vector<int>(1, INT_LIT(res));

}

// Set initial phase for a given variable
void
Maplesat::setPhase(const int var, const bool phase)
{
    solver->setPolarity(var - 1, phase ? true : false);
}

// Bump activity for a given variable
void
Maplesat::bumpVariableActivity(const int var, const int times)
{
//    for(int i = 0; i < times; i++) {
//        solver->varBumpActivity(var - 1, 1);
//        //TODO: work only for VSIDS
//    }
}

// Interrupt the SAT solving, so it can be started again with new assumptions
void
Maplesat::setSolverInterrupt()
{
    stopSolver = true;

    solver->interrupt();
}

void
Maplesat::unsetSolverInterrupt()
{
    stopSolver = false;

    solver->clearInterrupt();
}

// Diversify the solver
void
Maplesat::diversify(int id)
{
    solver->random_seed = (double)id;
//    if (id % 2) {
//        solver->VSIDS = true;
//    } else {
//        solver->VSIDS = false;
//    }
}

// Solve the formula with a given set of assumptions
// return 10 for SAT, 20 for UNSAT, 0 for UNKNOWN
SatResult
Maplesat::solve(const vector<int> & cube, int nConflicts)
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

    if (res == l_False) {
        if ( solver->conflict.size() > 0 ) {
            ClauseExchange *finalConf = getFinalAnalysis();
            clausesToExport.addClause(finalConf);
        }
        return UNSAT;
    }

    return UNKNOWN;
}

void
Maplesat::addClause(ClauseExchange * clause)
{
    clausesToAdd.addClause(clause);

    setSolverInterrupt();
}

void
Maplesat::addLearnedClause(ClauseExchange * clause)
{
    if (clause->size == 1) {
        unitsToImport.addClause(clause);
    } else {
        clausesToImport.addClause(clause);
    }
}

void
Maplesat::addClauses(const vector<ClauseExchange *> & clauses)
{
    clausesToAdd.addClauses(clauses);

    setSolverInterrupt();
}

void
Maplesat::addInitialClauses(const vector<ClauseExchange *> & clauses)
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
Maplesat::addLearnedClauses(const vector<ClauseExchange *> & clauses)
{
    for (size_t i = 0; i < clauses.size(); i++) {
        addLearnedClause(clauses[i]);
    }
}

void
Maplesat::getLearnedClauses(vector<ClauseExchange *> & clauses)
{
    clausesToExport.getClauses(clauses);
}

void
Maplesat::increaseClauseProduction()
{
    lbdLimit++;
}

void
Maplesat::decreaseClauseProduction()
{
    if (lbdLimit > 2) {
        lbdLimit--;
    }
}

SolvingStatistics
Maplesat::getStatistics()
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
Maplesat::getModel()
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
Maplesat::getFinalAnalysis()
{
   ClauseExchange * out_cls =
      ClauseManager::allocClause(solver->conflict.size());

   out_cls->lbd = 2;
   out_cls->from = id;

   for (int i = 0; i < solver->conflict.size(); i++) {
      out_cls->lits[i] = INT_LIT(solver->conflict[i]);
   }

   return out_cls;
}
