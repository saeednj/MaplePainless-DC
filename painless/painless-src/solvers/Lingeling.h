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

#include "../clauses/ClauseBuffer.h"
#include "../solvers/SolverInterface.h"
#include "../utils/SatUtils.h"
#include "../utils/Threading.h"

#include <atomic>

using namespace std;

struct LGL;

/// Instance of a Lingeling solver
class Lingeling: public SolverInterface
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
   Lingeling(int id);
   
   /// Copy constructor.
   Lingeling(const Lingeling & other, int id);

   /// Destructor.
   ~Lingeling();

protected:
   /// Pointer to a Lingeling solver
   LGL * solver;

   /// Boolean used to determine if the resolution should stop 
   atomic<int> stopSolver;

   /// LBD limit used to share clauses.
   atomic<int> glueLimit;
   
   /// Buffer used to add permanent clauses.
   ClauseBuffer clausesToAdd;

   /// Buffer used to import units.
   ClauseBuffer unitsToImport;

   /// Buffer used to import clauses.
   ClauseBuffer clausesToImport;

   /// Buffer used to export clauses (units included).
   ClauseBuffer clausesToExport;

   /// Size of the unit array used by Lingeling.
   size_t unitsBufferSize;
   
   /// Size of the clauses array used by Lingeling.
   size_t clsBufferSize;

   /// Unit array used by Lingeling.
   int * unitsBuffer;

   /// Clauses array used by Lingeling.
   int * clsBuffer;
   
   /// Termination callback.
   friend int  termCallback(void* solverPtr);

   /// Callback to export units.
   friend void produceUnit (void* sp, int lit);

   /// Callback to export clauses.
   friend void produce     (void* sp, int* cls, int glue);

   /// Callback to import units.
   friend void consumeUnits(void* sp, int** start, int** end);

   /// Callback to import clauses.
   friend void consumeCls  (void* sp, int** clause, int* glue);
};
