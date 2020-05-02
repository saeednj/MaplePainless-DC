#pragma once

#include "../clauses/ClauseBuffer.h"
#include "../solvers/SolverInterface.h"
#include "../utils/Threading.h"

using namespace std;

// Some forward declatarations for MapleSAT
namespace MapleSAT
{
	class SimpSolver;
	class Lit;
	template<class T> class vec;
}

/// Instance of a MapleSAT solver
class Maplesat : public SolverInterface
{
public:
   /// Load formula from a given dimacs file, return false if failed.
   bool loadFormula(const char* filename);

   /// Get the number of variables of the current resolution.
   int getVariablesCount();

   /// Get a variable suitable for search splitting.
   vector<int> getDivisionVariables(int k = 1);

   /// Set initial phase for a given variable.
   void setPhase(const int var, const bool phase);

   /// Bump activity of a given variable.
   void bumpVariableActivity(const int var, const int times);

   /// Interrupt resolution, solving cannot continue until interrupt is unset.
   void setSolverInterrupt();
   
   /// Remove the SAT solving interrupt request.
   void unsetSolverInterrupt();

   /// Solve the formula with a given cube.
   SatResult solve(const vector<int> & cube, int nConflicts = -1);

   /// Add a permanent clause to the formula.
   void addClause(ClauseExchange * clause);
   
   /// Add a list of permanent clauses to the formula.
   void addClauses(const vector<ClauseExchange *> & clauses);
   
   /// Add a list of initial clauses to the formula.
   void addInitialClauses(const vector<ClauseExchange *> & clauses);

   /// Add a learned clause to the formula.
   void addLearnedClause(ClauseExchange * clause);
   
   /// Add a list of learned clauses to the formula.
   void addLearnedClauses(const vector<ClauseExchange *> & clauses);

   /// Get a list of learned clauses.
   void getLearnedClauses(vector<ClauseExchange *> & clauses);

   /// Request the solver to produce more clauses.
   void increaseClauseProduction();
   
   /// Request the solver to produce less clauses.
   void decreaseClauseProduction();

   /// Get solver statistics.
   SolvingStatistics getStatistics();
   
   /// Return the model in case of SAT result.
   vector<int> getModel();
   
   /// Return the clause computed by the final conflict analysis.
   ClauseExchange * getFinalAnalysis();

   /// Native diversification.
   void diversify(int id);

   /// Constructor.
   Maplesat(int id);

   /// Copy Constructor
   Maplesat(const Maplesat & other, int id);

   /// Destructor.
   virtual ~Maplesat();

protected:
   /// Pointer to a Maplesat solver.
   MapleSAT::SimpSolver * solver;

   /// Buffer used to import units.
   ClauseBuffer unitsToImport;

   /// Buffer used to import clauses.
   ClauseBuffer clausesToImport;

   /// Buffer used to export clauses (units included).
   ClauseBuffer clausesToExport;

   /// Buffer used to add permanent clauses.
   ClauseBuffer clausesToAdd;
   
   /// Size limit used to share clauses.
   atomic<int> lbdLimit;
   
   /// Used to stop or continue the resolution.
   atomic<bool> stopSolver;
   
   /// Callback to import units.
   friend MapleSAT::Lit  maplesatImportUnit(void *);

   /// Callback to import clauses.
   friend bool maplesatImportClause(void *, int *, MapleSAT::vec<MapleSAT::Lit> &);

   /// Callback to export clauses (units included).
   friend void maplesatExportClause(void *, int, MapleSAT::vec<MapleSAT::Lit> &);
};
