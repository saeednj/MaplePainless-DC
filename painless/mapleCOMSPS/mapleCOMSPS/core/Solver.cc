/***************************************************************************************[Solver.cc]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

Chanseok Oh's MiniSat Patch Series -- Copyright (c) 2015, Chanseok Oh

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <math.h>
#include <signal.h>
#include <unistd.h>

#include "mapleCOMSPS/mtl/Sort.h"
#include "mapleCOMSPS/core/Solver.h"

#include "mapleCOMSPS/core/classifier.h"

#include <algorithm>
#include <unordered_set>

using namespace MapleCOMSPS;

#ifdef BIN_DRUP
int Solver::buf_len = 0;
unsigned char Solver::drup_buf[2 * 1024 * 1024];
unsigned char* Solver::buf_ptr = drup_buf;
#endif

//=================================================================================================
// Options:


static const char* _cat = "CORE";

static DoubleOption  opt_step_size         (_cat, "step-size",   "Initial step size",                             0.40,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_step_size_dec     (_cat, "step-size-dec","Step size decrement",                          0.000001, DoubleRange(0, false, 1, false));
static DoubleOption  opt_min_step_size     (_cat, "min-step-size","Minimal step size",                            0.06,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.80,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_clause_decay      (_cat, "cla-decay",   "The clause activity decay factor",              0.999,    DoubleRange(0, false, 1, false));
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq",    "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static DoubleOption  opt_random_seed       (_cat, "rnd-seed",    "Used by the random variable selection",         91648253, DoubleRange(0, false, HUGE_VAL, false));
static IntOption     opt_ccmin_mode        (_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));
static IntOption     opt_phase_saving      (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init",    "Randomize the initial activity", false);
static IntOption     opt_restart_first     (_cat, "rfirst",      "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption  opt_restart_inc       (_cat, "rinc",        "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));
static DoubleOption  opt_garbage_frac      (_cat, "gc-frac",     "The fraction of wasted memory allowed before a garbage collection is triggered",  0.20, DoubleRange(0, false, HUGE_VAL, false));

static BoolOption    opt_jw_pol      (_cat, "jw-pol",    "Jeroslow-Wang based polarity initialization", false);
static BoolOption    opt_jw_act      (_cat, "jw-act",    "Jeroslow-Wang based activity initialization", false);

//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :
    // Parallel parameters:
    //
    issuer                 (NULL)
  , exportClauseCallback   (NULL)
  , importUnitCallback     (NULL)
  , importClauseCallback   (NULL)

    // Parameters (user settable):
    //
  , drup_file        (NULL)
  , verbosity        (0)
  , step_size        (opt_step_size)
  , step_size_dec    (opt_step_size_dec)
  , min_step_size    (opt_min_step_size)
  , timer            (5000)
  , var_decay        (opt_var_decay)
  , clause_decay     (opt_clause_decay)
  , random_var_freq  (opt_random_var_freq)
  , random_seed      (opt_random_seed)
  , VSIDS            (false)
  , ccmin_mode       (opt_ccmin_mode)
  , phase_saving     (opt_phase_saving)
  , rnd_pol          (false)
  , rnd_init_act     (opt_rnd_init_act)
  , garbage_frac     (opt_garbage_frac)
  , restart_first    (opt_restart_first)
  , restart_inc      (opt_restart_inc)

    // Parameters (the rest):
    //
  , learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

    // Parameters (experimental):
    //
  , learntsize_adjust_start_confl (100)
  , learntsize_adjust_inc         (1.5)

    // Statistics: (formerly in 'SolverStats')
    //
  , solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0), conflicts_VSIDS(0)
  , dec_vars(0), clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)

  , ok                 (true)
  , cla_inc            (1)
  , var_inc            (1)
  , watches_bin        (WatcherDeleted(ca))
  , watches            (WatcherDeleted(ca))
  , qhead              (0)
  , simpDB_assigns     (-1)
  , simpDB_props       (0)
  , order_heap_CHB     (VarOrderLt(activity_CHB))
  , order_heap_VSIDS   (VarOrderLt(activity_VSIDS))
  , progress_estimate  (0)
  , remove_satisfied   (true)
    ,usePR(false)
    ,useFlip(false)

  , core_lbd_cut       (3)
  , tier2_lbd_cut      (6)
  , global_lbd_sum     (0)
  , lbd_queue          (50)
  , next_T2_reduce     (10000)
  , next_L_reduce      (15000)
  
  , counter            (0)

    // Resource constraints:
    //
  , conflict_budget    (-1)
  , propagation_budget (-1)
  , asynch_interrupt   (false)
  , jw_pol (opt_jw_pol)
  , jw_act (opt_jw_act)
  , strucFeatureFlag(0)
  , BMM(0)


{}

Solver::Solver(const Solver &s) :
  // Parallel parameters:
    //
    issuer                 (s.issuer)
  , exportClauseCallback   (s.exportClauseCallback)
  , importUnitCallback     (s.importUnitCallback)
  , importClauseCallback   (s.importClauseCallback)

    // Parameters (user settable):
    //
  , drup_file        (s.drup_file)
  , verbosity        (s.verbosity)
  , step_size        (s.step_size)
  , step_size_dec    (s.step_size_dec)
  , min_step_size    (s.min_step_size)
  , timer            (s.timer)
  , var_decay        (s.var_decay)
  , clause_decay     (s.clause_decay)
  , random_var_freq  (s.random_var_freq)
  , random_seed      (s.random_seed)
  , VSIDS            (s.VSIDS)
  , ccmin_mode       (s.ccmin_mode)
  , phase_saving     (s.phase_saving)
  , rnd_pol          (s.rnd_pol)
  , rnd_init_act     (s.rnd_init_act)
  , garbage_frac     (s.garbage_frac)
  , restart_first    (s.restart_first)
  , restart_inc      (s.restart_inc)
  , learntsize_factor(s.learntsize_factor)
  , learntsize_inc   (s.learntsize_inc)

  , learntsize_adjust_start_confl(s.learntsize_adjust_start_confl)
  , learntsize_adjust_inc(s.learntsize_adjust_inc)

  // Statistics: (formerly in 'SolverStats')
  //
  , solves(s.solves), starts(s.starts), decisions(s.decisions), rnd_decisions(s.rnd_decisions)
    , propagations(s.propagations), conflicts(s.conflicts), conflicts_VSIDS (s.conflicts_VSIDS)
  , dec_vars(s.dec_vars), clauses_literals(s.clauses_literals)
  , learnts_literals(s.learnts_literals), max_literals(s.max_literals), tot_literals(s.tot_literals)

  , ok(true)
  , cla_inc(s.cla_inc)
  , var_inc(s.var_inc)
  , watches(WatcherDeleted(ca))
  , watches_bin(WatcherDeleted(ca))
  , qhead(s.qhead)
  , simpDB_assigns(s.simpDB_assigns)
  , simpDB_props(s.simpDB_props)
  , order_heap_CHB(VarOrderLt(activity_CHB))
  , order_heap_VSIDS(VarOrderLt(activity_VSIDS))
  , progress_estimate(s.progress_estimate)
  , remove_satisfied(s.remove_satisfied)
  ,usePR(s.usePR)
  ,useFlip(s.useFlip)
  ,flipActivity(s.flipActivity)
  ,nbPropagations(s.nbPropagations)
  ,nbDecisionVar(s.nbDecisionVar)

  , core_lbd_cut       (s.core_lbd_cut)
  , tier2_lbd_cut      (s.tier2_lbd_cut)
  , global_lbd_sum     (s.global_lbd_sum)
  , lbd_queue          (s.lbd_queue)
  , next_T2_reduce     (s.next_T2_reduce)
  , next_L_reduce      (s.next_L_reduce)
  
  , counter            (s.counter)
  
  // Resource constraints:
    //
  , conflict_budget    (s.conflict_budget)
  , propagation_budget (s.propagation_budget)
  , asynch_interrupt   (s.asynch_interrupt)
  , jw_pol (s.jw_pol)
  , jw_act (s.jw_act)
  , strucFeatureFlag(s.strucFeatureFlag)
  , BMM(s.BMM)


{
  // Copy clauses.
  s.ca.copyTo(ca);
  ca.extra_clause_field = s.ca.extra_clause_field;
  
  s.binary.copyTo(binary);
  s.ternary.copyTo(ternary);
  s.numLearnt.copyTo(numLearnt);
  s.numAssigned.copyTo(numAssigned);
  
      // Copy all search vectors
   s.watches.copyTo(watches);
   s.watches_bin.copyTo(watches_bin);
   s.assigns.memCopyTo(assigns);
   s.vardata.memCopyTo(vardata);
   s.activity_CHB.memCopyTo(activity_CHB);
   s.activity_VSIDS.memCopyTo(activity_VSIDS);
   s.activity_lits.memCopyTo(activity_lits);
   s.seen.memCopyTo(seen);
   s.analyze_stack.memCopyTo(analyze_stack);
   s.analyze_toclear.memCopyTo(analyze_toclear);
   s.seen2.memCopyTo(seen2);
   s.add_tmp.memCopyTo(add_tmp);
   s.add_oc.memCopyTo(add_oc);
   s.polarity.memCopyTo(polarity);
   s.decision.memCopyTo(decision);
   s.trail.memCopyTo(trail);
   s.order_heap_CHB.copyTo(order_heap_CHB);
   s.order_heap_VSIDS.copyTo(order_heap_VSIDS);
   s.clauses.memCopyTo(clauses);
   s.learnts_core.memCopyTo(learnts_core);
   s.learnts_tier2.memCopyTo(learnts_tier2);
   s.learnts_local.memCopyTo(learnts_local);
   s.picked.memCopyTo(picked);
   s.conflicted.memCopyTo(conflicted);
   s.almost_conflicted.memCopyTo(almost_conflicted);
   s.parameters.memCopyTo(parameters);
   s.updatedParams.memCopyTo(updatedParams);
#ifdef ANTI_EXPLORATION
   s.canceled.memCopyTo(canceled);
#endif

}

Solver::~Solver()
{
}


//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar)
{
    int v = nVars();
    watches_bin.init(mkLit(v, false));
    watches_bin.init(mkLit(v, true ));
    watches  .init(mkLit(v, false));
    watches  .init(mkLit(v, true ));
    assigns  .push(l_Undef);
    vardata  .push(mkVarData(CRef_Undef, 0));
    activity_CHB  .push(0);
    activity_VSIDS.push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    activity_lits.push(0);
    activity_lits.push(0);
    if(useFlip) {
       flipActivity.push_back(0);
    }
    if(usePR){
      nbPropagations.push_back(0);
      nbDecisionVar.push_back(0);
    }
    picked.push(0);
    conflicted.push(0);
    almost_conflicted.push(0);
#ifdef ANTI_EXPLORATION
    canceled.push(0);
#endif

    seen     .push(0);
    seen2    .push(0);
    polarity .push(sign);
    decision .push();
    trail    .capacity(v+1);
    setDecisionVar(v, dvar);
    binary.push(0);
    ternary.push(0);
    numLearnt.push(0);
    numAssigned.push(0);
    return v;
}


bool Solver::addClause_(vec<Lit>& ps)
{
    assert(decisionLevel() == 0);
    if (!ok) return false;

    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);
    Lit p; int i, j;

    if (drup_file){
        add_oc.clear();
        for (int i = 0; i < ps.size(); i++) add_oc.push(ps[i]); }

    for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
        if (value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if (value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    ps.shrink(i - j);

    if (drup_file && i != j){
#ifdef BIN_DRUP
        binDRUP('a', ps, drup_file);
        binDRUP('d', add_oc, drup_file);
#else
        for (int i = 0; i < ps.size(); i++)
            fprintf(drup_file, "%i ", (var(ps[i]) + 1) * (-2 * sign(ps[i]) + 1));
        fprintf(drup_file, "0\n");

        fprintf(drup_file, "d ");
        for (int i = 0; i < add_oc.size(); i++)
            fprintf(drup_file, "%i ", (var(add_oc[i]) + 1) * (-2 * sign(add_oc[i]) + 1));
        fprintf(drup_file, "0\n");
#endif
    }

    if (ps.size() == 0)
        return ok = false;
    else if (ps.size() == 1){
        uncheckedEnqueue(ps[0]);
        return ok = (propagate() == CRef_Undef);
    }else{
        CRef cr = ca.alloc(ps, false);
        clauses.push(cr);
        attachClause(cr);
    }

    return true;
}

vector< map<int, double> > Solver::createVIG() {
   vector< map<int, double> > graph;
   graph.resize(nVars(), map<int, double>());

   for(int i = 0; i < clauses.size(); i++) {
      Clause & c = ca[clauses[i]];
      double weight = 1.0;
      double comb = c.size() * (c.size() - 1) / 2;
      weight /= comb;

      for(int j = 0; j < c.size(); j++) {
         int var0 = var(c[j]);

         for(int k = j + 1; k < c.size(); k++) {
            int var1 = var(c[k]);
            if(var0 != var1) {
               graph[var0][var1] += weight;
               graph[var1][var0] += weight;
            }
         }
      }
   }

   //for(int i = 0; i < learnts_core.size(); i++) {
   //   Clause & c = ca[learnts_core[i]];
   //   double weight = 1.0;
   //   double comb = c.size() * (c.size() - 1) / 2;
   //   weight /= comb;

   //   for(int j = 0; j < c.size(); j++) {
   //      int var0 = var(c[j]);

   //      for(int k = j + 1; k < c.size(); k++) {
   //         int var1 = var(c[k]);
   //         if(var0 != var1) {
   //            graph[var0][var1] += weight;
   //            graph[var1][var0] += weight;
   //         }
   //      }
   //   }
   //}


   //for(int i = 0; i < learnts_tier2.size(); i++) {
   //   Clause & c = ca[learnts_tier2[i]];
   //   double weight = 1.0;
   //   double comb = c.size() * (c.size() - 1) / 2;
   //   weight /= comb;

   //   for(int j = 0; j < c.size(); j++) {
   //      int var0 = var(c[j]);

   //      for(int k = j + 1; k < c.size(); k++) {
   //         int var1 = var(c[k]);
   //         if(var0 != var1) {
   //            graph[var0][var1] += weight;
   //            graph[var1][var0] += weight;
   //         }
   //      }
   //   }
   //}

   //for(int i = 0; i < learnts_local.size(); i++) {
   //   Clause & c = ca[learnts_local[i]];
   //   double weight = 1.0;
   //   double comb = c.size() * (c.size() - 1) / 2;
   //   weight /= comb;

   //   for(int j = 0; j < c.size(); j++) {
   //      int var0 = var(c[j]);

   //      for(int k = j + 1; k < c.size(); k++) {
   //         int var1 = var(c[k]);
   //         if(var0 != var1) {
   //            graph[var0][var1] += weight;
   //            graph[var1][var0] += weight;
   //         }
   //      }
   //   }
   //}

   return graph;
}

void Solver::importUnitClauses() {
   assert(decisionLevel() == 0);

   if (importUnitCallback == NULL)
      return;

   Lit l;
   while ((l = importUnitCallback(issuer)) != lit_Undef) {
      if (value(var(l)) == l_Undef) {
         uncheckedEnqueue(l);
      }
   }
}

bool Solver::importClauses() {
    assert(decisionLevel() == 0);

    int lbd, k, l;
    bool alreadySat;
    while (importClauseCallback(issuer, &lbd, importedClause)) {
        alreadySat = false;
        // Simplify clause before add
        for (k = l = 0; k < importedClause.size(); k++) {
            if (value(importedClause[k]) == l_True) {
                alreadySat = true;
                break;
            } else if (value(importedClause[k]) == l_Undef) {
                importedClause[l++] = importedClause[k];
            }
        }
        importedClause.shrink(k - l);

        if (alreadySat) {
            importedClause.clear();
            continue;
        }

        if (importedClause.size() == 0) {
           return false;
        } else if (importedClause.size() == 1) {
            uncheckedEnqueue(importedClause[0]);
        } else {
            CRef cr = ca.alloc(importedClause, true);
            if (importedClause.size() == 2)
                lbd = lbd > 2 ? 2 : lbd;
            ca[cr].set_lbd(lbd);
            if (lbd <= core_lbd_cut) {
                learnts_core.push(cr);
                ca[cr].mark(CORE);
            } else if (lbd <= 6) {
                learnts_tier2.push(cr);
                ca[cr].mark(TIER2);
                ca[cr].touched() = conflicts;
            } else {
                learnts_local.push(cr);
                claBumpActivity(ca[cr]);
            }
            attachClause(cr);
        }
        importedClause.clear();
    }

    return true;
}
//void Solver::importClauses() {
//   assert(decisionLevel() == 0);
//
//   if (importClauseCallback == NULL)
//      return;
//
//   int lbd;
//   while (importClauseCallback(issuer, &lbd, importedClause)) {
//      CRef cr = ca.alloc(importedClause, true);
//      ca[cr].set_lbd(lbd);
//      if (lbd <= core_lbd_cut) {
//         learnts_core.push(cr);
//         ca[cr].mark(CORE);
//      } else if (lbd <= tier2_lbd_cut) {
//         learnts_tier2.push(cr);
//         ca[cr].mark(TIER2);
//         ca[cr].touched() = conflicts;
//      } else {
//         learnts_local.push(cr);
//         claBumpActivity(ca[cr]);
//      }
//      attachClause(cr);
//      importedClause.clear();
//   }
//}


void Solver::attachClause(CRef cr) {
    const Clause& c = ca[cr];
    assert(c.size() > 1);
    OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = c.size() == 2 ? watches_bin : watches;
    ws[~c[0]].push(Watcher(cr, c[1]));
    ws[~c[1]].push(Watcher(cr, c[0]));
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size(); }


void Solver::detachClause(CRef cr, bool strict) {
    const Clause& c = ca[cr];
    assert(c.size() > 1);
    OccLists<Lit, vec<Watcher>, WatcherDeleted>& ws = c.size() == 2 ? watches_bin : watches;
    
    if (strict){
        remove(ws[~c[0]], Watcher(cr, c[1]));
        remove(ws[~c[1]], Watcher(cr, c[0]));
    }else{
        // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
        ws.smudge(~c[0]);
        ws.smudge(~c[1]);
    }

    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size(); }


void Solver::removeClause(CRef cr) {
    Clause& c = ca[cr];

    if (drup_file){
        if (c.mark() != 1){
#ifdef BIN_DRUP
            binDRUP('d', c, drup_file);
#else
            fprintf(drup_file, "d ");
            for (int i = 0; i < c.size(); i++)
                fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
            fprintf(drup_file, "0\n");
#endif
        }else
            printf("c Bug. I don't expect this to happen.\n");
    }

    detachClause(cr);
    // Don't leave pointers to free'd memory!
    if (locked(c)){
        Lit implied = c.size() != 2 ? c[0] : (value(c[0]) == l_True ? c[0] : c[1]);
        vardata[var(implied)].reason = CRef_Undef; }
    c.mark(1); 
    ca.free(cr);
}


bool Solver::satisfied(const Clause& c) const {
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; }


// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
    if (decisionLevel() > level){
        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
            Var      x  = var(trail[c]);

            if (!VSIDS){
                uint32_t age = conflicts - picked[x];
                if (age > 0){
                    double adjusted_reward = ((double) (conflicted[x] + almost_conflicted[x])) / ((double) age);
                    double old_activity = activity_CHB[x];
                    activity_CHB[x] = step_size * adjusted_reward + ((1 - step_size) * old_activity);

                    int idx = toInt(trail[c]);
                    double old_lit_activity = activity_lits[idx];
                    activity_lits[idx] = step_size * adjusted_reward + ((1 - step_size) * old_lit_activity);

                    if (order_heap_CHB.inHeap(x)){
                        if (activity_CHB[x] > old_activity)
                            order_heap_CHB.decrease(x);
                        else
                            order_heap_CHB.increase(x);
                    }
                }
#ifdef ANTI_EXPLORATION
                canceled[x] = conflicts;
#endif
            }

            assigns [x] = l_Undef;
            if (phase_saving > 1 || (phase_saving == 1) && c > trail_lim.last()) {
                bool b = sign(trail[c]);
                if(useFlip && b!=polarity[x]) flipActivity[x]++;
                polarity[x] = b;
            }
            insertVarOrder(x); }
        qhead = trail_lim[level];
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);
    } }


//=================================================================================================
// Major methods:

// activity
// LRBproduct
// flip
// avg_propagation
// num_decided
// num_assigned
// num_learnt
// decision_level
// num_binary_clause
// num_ternary_clause
vector<int> Solver::pickSplittingVariables(int k)
{
    vector<int> r;

    unordered_set<int> assumSet, pickedSet;
    for( int i=0; i<assumptions.size(); i++ )
        assumSet.insert(var(assumptions[i]));


    float features[2*NUM_FEATURES];
    if (!strucFeatureFlag)
    {
        strucFeatureFlag = 1;
        for( int i=0; i<nVars(); i++ )
            binary[i] = ternary[i] = 0;
        for( int i=0; i<clauses.size(); i++ )
        {
            CRef cr = clauses[i];
            Clause& c = ca[cr];
            if ( c.size() == 2 && value(c[0]) != l_True && value(c[1]) != l_True )
            {
                binary[var(c[0])]++;
                binary[var(c[1])]++;
            }
            if ( c.size() == 3 && value(c[0]) != l_True && value(c[1]) != l_True && value(c[2]) != l_True )
            {
                ternary[var(c[0])]++;
                ternary[var(c[1])]++;
                ternary[var(c[2])]++;
            }
        }
    }
    int first = 0;
    for( int i=0; i<k; i++ )
    {
        while( (assumSet.count(first) > 0 || pickedSet.count(first) > 0) && first < nVars() )
            first++;
        if ( first == nVars() ) return r;
        int best = first;
        getFeatures(features, NUM_FEATURES, best);
        for( int v=first+1; v<nVars(); v++ )
        {
            if ( assumSet.count(v) > 0 ) continue;

            getFeatures(features, 0, v);
            if ( predict(features) ) // is 'v' better than 'best' ?
            {
                best = v;
                for( int j=0; j<NUM_FEATURES; j++ )
                    features[j+NUM_FEATURES] = features[j];
            }
        }
        pickedSet.insert(best);
        r.push_back(best);
    }
    return r;
}

void Solver::getFeatures(float* features, int offset, int v)
{
    features[offset  ] = activity_VSIDS[v];
    double n_val = activity_lits[2*v];
    double p_val = activity_lits[2*v+1];
    double lit_value = n_val * p_val * 1024 + n_val + p_val;
    features[offset+1] = lit_value;
    features[offset+2] = flipActivity[v];
    features[offset+3] = nbDecisionVar[v] == 0 ? 0 : nbPropagations[v] * 1.0 / nbDecisionVar[v];
    features[offset+4] = nbDecisionVar[v];
    features[offset+5] = numAssigned[v];
    features[offset+5] = numLearnt[v];
    features[offset+6] = level(v);
    features[offset+7] = binary[v];
    features[offset+8] = ternary[v];
}

// Returns 1 if first is better
//     and 0 if second is better
//int Solver::splittingPredicate(float* features, int first, int second)
//{
//    return predict();
//}



Lit Solver::pickBranchLit()
{
    Var next = var_Undef;
    Heap<VarOrderLt>& order_heap = VSIDS ? order_heap_VSIDS : order_heap_CHB;

    // Random decision:
    /*if (drand(random_seed) < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(random_seed,order_heap.size())];
        if (value(next) == l_Undef && decision[next])
            rnd_decisions++; }*/

    // Activity based decision:
    while (next == var_Undef || value(next) != l_Undef || !decision[next])
        if (order_heap.empty())
            return lit_Undef;
        else{
#ifdef ANTI_EXPLORATION
            if (!VSIDS){
                Var v = order_heap_CHB[0];
                uint32_t age = conflicts - canceled[v];
                while (age > 0){
                    double decay = pow(0.95, age);
                    activity_CHB[v] *= decay;
                    activity_lits[2*v] *= decay;
                    activity_lits[2*v+1] *= decay;
                    if (order_heap_CHB.inHeap(v))
                        order_heap_CHB.increase(v);
                    canceled[v] = conflicts;
                    v = order_heap_CHB[0];
                    age = conflicts - canceled[v];
                }
            }
#endif
            next = order_heap.removeMin();
        }

    return mkLit(next, polarity[next]);
}

void Solver::pickBranchLit(int n, std::vector<Lit>& lits){
  /*int cpt=0;
  (VSIDS?order_heap_VSIDS:order_heap_CHB).copyTo(order_heap_copy);
  Var next=var_Undef;

  while(cpt++<n && !order_heap_copy.empty()){
    next = order_heap_copy.removeMin();
    if(next!=var_Undef && value(next)==l_Undef && decision[next])
      lits.push_back(mkLit(next, polarity[next]));
      }*/
  printf("Fatal: multi pickBranchLit not implemented yet\n");
  abort();
}

vector<int> Solver::pickSplittingVariablesUsingLRBLit(int k)
{
    Var next = 0;
    double max_val = -1;

    for(Var v = 0; v < nVars(); v++) {
        if (decision[v] && value(v) == l_Undef) {
            double n_val = activity_lits[2*v];
            double p_val = activity_lits[2*v+1];
            double lit_value = n_val * p_val * 1024 + n_val + p_val;
            if (lit_value > max_val) {
                max_val = lit_value;
                next = v;
            }
        }
    }

    return vector<int>(1, max_val == -1 ? var_Undef : next);
}


Lit Solver::pickBranchLitUsingBayesian() {
    Var next = 0;
    double min_diff = 1000000000;

    //printf("DBG: learnts_core.size: %d\n", learnts_core.size());
    //fflush(stdout);
    for( int i=0; i<learnts_core.size(); i++ )
    {
        Clause& c = ca[learnts_core[i]];
        if ( c.size() <= 3 )
            bayesian_update(c);
    }
    int max_flip = -1;
    for( Var v=0; v<nVars(); v++ ) {
        double a = parameters[v].a;
        double b = parameters[v].b;
        if ( a < 3*b || b < 3*a )
        {
            if ( decision[v] && value(v) == l_Undef && flipActivity[v] > max_flip )
            {
                max_flip = flipActivity[v];
                next = v;
            }
        }
        //double diff = fabs(parameters[v].a - parameters[v].b);
        //if ( decision[v] && value(v) == l_Undef && diff < min_diff ) {
        //    next = v;
        //    min_diff = diff;
        //}
    }
    if ( max_flip == -1 )
    {
        for( Var v=0; v<nVars(); v++ ) {
            if ( decision[v] && value(v) == l_Undef && flipActivity[v] > max_flip )
            {
                max_flip = flipActivity[v];
                next = v;
            }
        }
    }
    return mkLit(next, polarity[next]);
}

Lit Solver::pickBranchLitUsingFlipActivity() {
  assert(useFlip);
  Var next = 0;
  int max_flips = -1;

   for(Var v = 0; v < nVars(); v++) {
      if(flipActivity[v] > max_flips && decision[v] && value(v) == l_Undef) {
         max_flips = flipActivity[v];
         next      = v;
      }
   }

   return max_flips == -1 ? lit_Undef : mkLit(next, polarity[next]);
}

void Solver::pickBranchLitUsingFlipActivity(int n, std::vector<Lit>& lits){
  assert(useFlip);
  int cpt=0;
  Var next;
  int val;
  std::vector<int> copy=std::vector<int>(flipActivity);

  do{
    next=0;
    val=-1;
    for(size_t i=0;i<copy.size();i++){
      if(copy[i]>val && decision[i] && value(i) == l_Undef){
	val=copy[i];
	next=i;
      }
    }
    if(val==-1) break;
    lits.push_back(mkLit(next, rnd_pol ? drand(random_seed) < 0.5 : polarity[next]));
    copy[next]=-1;

  }while(cpt++<n);
}

Lit Solver::pickBranchLitUsingPropagationRate() {
  assert(usePR);
  Var next      = 0;
  double max_pr = -1.0;
  double tmp_pr;

   for(Var v = 0; v < nVars(); v++) {
      if (nbDecisionVar[v] == 0) {
         tmp_pr = 0.0;
      } else {
         tmp_pr = nbPropagations[v] * 1.0 / nbDecisionVar[v];
      }
      if(nbDecisionVar[v] != 0 && tmp_pr > max_pr && decision[v] && value(v) == l_Undef) {
         max_pr = tmp_pr;
         next   = v;
      }
   }

   fill(nbDecisionVar.begin(), nbDecisionVar.end(), 0);
   fill(nbPropagations.begin(), nbPropagations.end(), 0);

   return max_pr == -1.0 ? lit_Undef : mkLit(next, polarity[next]);
}

void Solver::pickBranchLitUsingPropagationRate(int n,std::vector<Lit>& lits){
  assert(usePR);
  int cpt=0;
  Var next;
  double val,tmp;
  std::vector<int> copy=std::vector<int>(nbDecisionVar);


  do{
    next=0;
    val=-1.0;
    for(size_t i=0;i<copy.size();i++){
      if(copy[i]!=0 && (tmp=(nbPropagations[i]*1.0/copy[i]))>val && decision[i] && value(i) == l_Undef){
	val=tmp;
	next=i;
      }
    }
    if(val==-1.0) break;
    lits.push_back(mkLit(next, rnd_pol ? drand(random_seed) < 0.5 : polarity[next]));
    copy[next]=0;

  }while(cpt++<n);
}



/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|  
|  Description:
|    Analyze conflict and produce a reason clause.
|  
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|  
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the 
|        rest of literals. There may be others from the same level though.
|  
|________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, vec<Lit>& out_learnt, int & out_btlevel,
                     int & out_lbd)
{
    int pathC    = 0;
    Lit p        = lit_Undef;

    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;

    do{
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied lit.
        if (p != lit_Undef && c.size() == 2 && value(c[0]) == l_False){
            assert(value(c[1]) == l_True);
            Lit tmp = c[0];
            c[0] = c[1], c[1] = tmp; }

        // Update LBD if improved.
        if (c.learnt() && c.mark() != CORE){
            int lbd = computeLBD(c);
            if (lbd < c.lbd()){
                if (c.lbd() <= 30) c.removable(false); // Protect once from reduction.
                c.set_lbd(lbd);
                if (lbd <= core_lbd_cut){
                    learnts_core.push(confl);
                    c.mark(CORE);
                }else if (lbd <= tier2_lbd_cut && c.mark() == LOCAL){
                    // Bug: 'cr' may already be in 'learnts_tier2', e.g., if 'cr' was demoted from TIER2
                    // to LOCAL previously and if that 'cr' is not cleaned from 'learnts_tier2' yet.
                    learnts_tier2.push(confl);
                    c.mark(TIER2); }
            }

            if (c.mark() == TIER2)
                c.touched() = conflicts;
            else if (c.mark() == LOCAL)
                claBumpActivity(c);
        }

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!seen[var(q)] && level(var(q)) > 0){
                if (VSIDS){
                    varBumpActivity(var(q), .5);
                    add_tmp.push(q);
                }else
                    conflicted[var(q)]++;
                seen[var(q)] = 1;
                if (level(var(q)) >= decisionLevel()){
                    pathC++;
                }else
                    out_learnt.push(q);
            }
        }
        
        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    }while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(analyze_toclear);
    if (ccmin_mode == 2){
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level))
                out_learnt[j++] = out_learnt[i];
        
    }else if (ccmin_mode == 1){
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[reason(var(out_learnt[i]))];
                for (int k = c.size() == 2 ? 0 : 1; k < c.size(); k++)
                    if (!seen[var(c[k])] && level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; }
            }
        }
    }else
        i = j = out_learnt.size();

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();

    for( int i=0; i<out_learnt.size(); i++ )
        numLearnt[var(out_learnt[i])]++;

    out_lbd = computeLBD(out_learnt);
    if (out_lbd <= tier2_lbd_cut && out_learnt.size() <= 30) // Try further minimization?
        if (binResMinimize(out_learnt))
            out_lbd = computeLBD(out_learnt); // Recompute LBD if minimized.

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level(var(p));
    }

    if (VSIDS){
        for (int i = 0; i < add_tmp.size(); i++){
            Var v = var(add_tmp[i]);
            if (level(v) >= out_btlevel - 1)
                varBumpActivity(v, 1);
        }
        add_tmp.clear();
    }else{
        seen[var(p)] = true;
        for(int i = out_learnt.size() - 1; i >= 0; i--){
            Var v = var(out_learnt[i]);
            CRef rea = reason(v);
            if (rea != CRef_Undef){
                const Clause& reaC = ca[rea];
                for (int i = 0; i < reaC.size(); i++){
                    Lit l = reaC[i];
                    if (!seen[var(l)]){
                        seen[var(l)] = true;
                        almost_conflicted[var(l)]++;
                        analyze_toclear.push(l); } } } } }

    for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
}


// Try further learnt clause minimization by means of binary clause resolution.
bool Solver::binResMinimize(vec<Lit>& out_learnt)
{
    // Preparation: remember which false variables we have in 'out_learnt'.
    counter++;
    for (int i = 1; i < out_learnt.size(); i++)
        seen2[var(out_learnt[i])] = counter;

    // Get the list of binary clauses containing 'out_learnt[0]'.
    const vec<Watcher>& ws = watches_bin[~out_learnt[0]];

    int to_remove = 0;
    for (int i = 0; i < ws.size(); i++){
        Lit the_other = ws[i].blocker;
        // Does 'the_other' appear negatively in 'out_learnt'?
        if (seen2[var(the_other)] == counter && value(the_other) == l_True){
            to_remove++;
            seen2[var(the_other)] = counter - 1; // Remember to remove this variable.
        }
    }

    // Shrink.
    if (to_remove > 0){
        int last = out_learnt.size() - 1;
        for (int i = 1; i < out_learnt.size() - to_remove; i++)
            if (seen2[var(out_learnt[i])] != counter)
                out_learnt[i--] = out_learnt[last--];
        out_learnt.shrink(to_remove);
    }
    return to_remove != 0;
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels)
{
    analyze_stack.clear(); analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(reason(var(analyze_stack.last())) != CRef_Undef);
        Clause& c = ca[reason(var(analyze_stack.last()))]; analyze_stack.pop();

        // Special handling for binary clauses like in 'analyze()'.
        if (c.size() == 2 && value(c[0]) == l_False){
            assert(value(c[1]) == l_True);
            Lit tmp = c[0];
            c[0] = c[1], c[1] = tmp; }

        for (int i = 1; i < c.size(); i++){
            Lit p  = c[i];
            if (!seen[var(p)] && level(var(p)) > 0){
                if (reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0){
                    seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }else{
                    for (int j = top; j < analyze_toclear.size(); j++)
                        seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }

    return true;
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|  
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict)
{
    out_conflict.clear();
    out_conflict.push(p);

    if (decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if (seen[x]){
            if (reason(x) == CRef_Undef){
                assert(level(x) > 0);
                out_conflict.push(~trail[i]);
            }else{
                Clause& c = ca[reason(x)];
                for (int j = c.size() == 2 ? 0 : 1; j < c.size(); j++)
                    if (level(var(c[j])) > 0)
                        seen[var(c[j])] = 1;
            }
            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
}


void Solver::uncheckedEnqueue(Lit p, CRef from)
{
    assert(value(p) == l_Undef);
    Var x = var(p);
    if (!VSIDS){
        picked[x] = conflicts;
        conflicted[x] = 0;
        almost_conflicted[x] = 0;
#ifdef ANTI_EXPLORATION
        uint32_t age = conflicts - canceled[var(p)];
        if (age > 0){
            double decay = pow(0.95, age);
            activity_CHB[var(p)] *= decay;
            activity_lits[2*var(p)] *= decay;
            activity_lits[1+2*var(p)] *= decay;
            if (order_heap_CHB.inHeap(var(p)))
                order_heap_CHB.increase(var(p));
        }
#endif
    }

    assigns[x] = lbool(!sign(p));
    numAssigned[x]++;
    vardata[x] = mkVarData(from, decisionLevel());
    trail.push_(p);
}


/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|  
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise CRef_Undef.
|  
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
CRef Solver::propagate()
{
    CRef    confl     = CRef_Undef;
    int     num_props = 0;
    int previousqhead = qhead;
    watches.cleanAll();
    watches_bin.cleanAll();
    Lit last = lit_Undef;
    if ( trail.size() > 0 )
    {
        last = trail[trail.size()-1];
    }

    while (qhead < trail.size()){
        Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
        vec<Watcher>&  ws  = watches[p];
        Watcher        *i, *j, *end;
        num_props++;

        vec<Watcher>& ws_bin = watches_bin[p];  // Propagate binary clauses first.
        for (int k = 0; k < ws_bin.size(); k++){
            Lit the_other = ws_bin[k].blocker;
            if (value(the_other) == l_False){
                confl = ws_bin[k].cref;
#ifdef LOOSE_PROP_STAT
                return confl;
#else
                goto ExitProp;
#endif
            }else if(value(the_other) == l_Undef)
                uncheckedEnqueue(the_other, ws_bin[k].cref);
        }

        for (i = j = (Watcher*)ws, end = i + ws.size();  i != end;){
            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (value(blocker) == l_True){
                *j++ = *i++; continue; }

            // Make sure the false literal is data[1]:
            CRef     cr        = i->cref;
            Clause&  c         = ca[cr];
            Lit      false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;
            assert(c[1] == false_lit);
            i++;

            // If 0th watch is true, then clause is already satisfied.
            Lit     first = c[0];
            Watcher w     = Watcher(cr, first);
            if (first != blocker && value(first) == l_True){
                *j++ = w; continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (value(c[k]) != l_False){
                    c[1] = c[k]; c[k] = false_lit;
                    watches[~c[1]].push(w);
                    goto NextClause; }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
            if (value(first) == l_False){
                confl = cr;
                qhead = trail.size();
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            }else
                uncheckedEnqueue(first, cr);

        NextClause:;
        }
        ws.shrink(i - j);
    }

    
 ExitProp:;
    if (usePR && num_props > 0)
    {
        nbPropagations[var(last)] += num_props;
        nbDecisionVar[var(last)]  += 1;
    }

    propagations += num_props;
    simpDB_props -= num_props;

    return confl;
}


/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|  
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt { 
    ClauseAllocator& ca;
    reduceDB_lt(ClauseAllocator& ca_) : ca(ca_) {}
    bool operator () (CRef x, CRef y) const { return ca[x].activity() < ca[y].activity(); }
};
void Solver::reduceDB()
{
    int     i, j;
    //if (local_learnts_dirty) cleanLearnts(learnts_local, LOCAL);
    //local_learnts_dirty = false;

    sort(learnts_local, reduceDB_lt(ca));

    int limit = learnts_local.size() / 2;
    for (i = j = 0; i < learnts_local.size(); i++){
        Clause& c = ca[learnts_local[i]];
        if (c.mark() == LOCAL)
            if (c.removable() && !locked(c) && i < limit)
                removeClause(learnts_local[i]);
            else{
                if (!c.removable()) limit++;
                c.removable(true);
                learnts_local[j++] = learnts_local[i]; }
    }
    learnts_local.shrink(i - j);

    checkGarbage();
}
void Solver::reduceDB_Tier2()
{
    int i, j;
    for (i = j = 0; i < learnts_tier2.size(); i++){
        Clause& c = ca[learnts_tier2[i]];
        if (c.mark() == TIER2)
            if (!locked(c) && c.touched() + 30000 < conflicts){
                learnts_local.push(learnts_tier2[i]);
                c.mark(LOCAL);
                //c.removable(true);
                c.activity() = 0;
                claBumpActivity(c);
            }else
                learnts_tier2[j++] = learnts_tier2[i];
    }
    learnts_tier2.shrink(i - j);
}


void Solver::removeSatisfied(vec<CRef>& cs)
{
    int i, j;
    for (i = j = 0; i < cs.size(); i++){
        Clause& c = ca[cs[i]];
        if (satisfied(c))
            removeClause(cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}

// TODO: very dirty and hackish.
void Solver::removeClauseHack(CRef cr, Lit watched0, Lit watched1)
{
    assert(ca[cr].size() >= 2);

    Clause& c = ca[cr];
    if (drup_file) // Hackish.
        if (c.mark() != 1){
#ifdef BIN_DRUP
            binDRUP('d', add_oc, drup_file); // 'add_oc' not 'c'.
#else
            for (int i = 0; i < add_oc.size(); i++)
                fprintf(drup_file, "%i ", (var(add_oc[i]) + 1) * (-2 * sign(add_oc[i]) + 1));
            fprintf(drup_file, "0\n");
#endif
        }else
            printf("c Bug: removeClauseHack(). I don't expect this to happen.\n");

    // TODO: dirty hack to exploit 'detachClause'. 'c' hasn't shrunk yet, so this will work fine.
    c[0] = watched0, c[1] = watched1;
    detachClause(cr);
    // Don't leave pointers to free'd memory!
    if (locked(c)){
        Lit implied = c.size() != 2 ? c[0] : (value(c[0]) == l_True ? c[0] : c[1]);
        vardata[var(implied)].reason = CRef_Undef; }
    c.mark(1);
    ca.free(cr);
}

// TODO: needs clean up.
void Solver::safeRemoveSatisfiedCompact(vec<CRef>& cs, unsigned valid_mark)
{
    int i, j, k, l;
    for (i = j = 0; i < cs.size(); i++){

        Clause& c = ca[cs[i]];
        if (c.mark() != valid_mark) continue;

        Lit c0 = c[0], c1 = c[1];
        if (drup_file){ // Remember the original clause before attempting to modify it.
            add_oc.clear();
            for (int i = 0; i < c.size(); i++) add_oc.push(c[i]); }

        // Remove false literals at the same time.
        for (k = l = 0; k < c.size(); k++)
            if (value(c[k]) == l_True){
                removeClauseHack(cs[i], c0, c1);
                goto NextClause; // Clause already satisfied; forget about it.
            }else if (value(c[k]) == l_Undef){
                c[l++] = c[k];
            }

        assert(1 < l && l <= k);

        // If became binary, we also need to migrate watchers. The easiest way is to allocate a new binary.
        if (l == 2 && k != 2){
            assert(add_tmp.size() == 0);
            add_tmp.push(c[0]); add_tmp.push(c[1]);
            bool learnt = c.learnt(); // Need a copy; see right below.
            int lbd = c.lbd();
            int m = c.mark();
            CRef cr = ca.alloc(add_tmp, learnt); // Caution! 'alloc' may invalidate the 'c' reference.
            if (learnt){
                if (m != CORE) learnts_core.push(cr);
                ca[cr].mark(CORE);
                ca[cr].set_lbd(lbd > 2 ? 2 : lbd); }
            attachClause(cr);

            if (drup_file){
#ifdef BIN_DRUP
                binDRUP('a', add_tmp, drup_file);
#else
                for (int i = 0; i < add_tmp.size(); i++)
                    fprintf(drup_file, "%i ", (var(add_tmp[i]) + 1) * (-2 * sign(add_tmp[i]) + 1));
                fprintf(drup_file, "0\n");
#endif
            }
            add_tmp.clear();

            removeClauseHack(cs[i], c0, c1);
            cs[j++] = cr; // Should be after 'removeClauseHack' because 'i' may be equal to 'j'.
            goto NextClause;
        }

        c.shrink(k - l); // FIXME: fix (statistical) memory leak.
        if (c.learnt()) learnts_literals -= (k - l);
        else            clauses_literals -= (k - l);

        if (drup_file && k != l){
#ifdef BIN_DRUP
            binDRUP('a', c, drup_file);
            binDRUP('d', add_oc, drup_file);
#else
            for (int i = 0; i < c.size(); i++)
                fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
            fprintf(drup_file, "0\n");

            fprintf(drup_file, "d ");
            for (int i = 0; i < add_oc.size(); i++)
                fprintf(drup_file, "%i ", (var(add_oc[i]) + 1) * (-2 * sign(add_oc[i]) + 1));
            fprintf(drup_file, "0\n");
#endif
        }

        cs[j++] = cs[i];
NextClause:;
    }
    cs.shrink(i - j);
}
//void Solver::safeRemoveSatisfied(vec<CRef>& cs, unsigned valid_mark)
//{
//    int i, j;
//    for (i = j = 0; i < cs.size(); i++){
//        Clause& c = ca[cs[i]];
//        if (c.mark() == valid_mark)
//            if (satisfied(c))
//                removeClause(cs[i]);
//            else
//                cs[j++] = cs[i];
//    }
//    cs.shrink(i - j);
//}

void Solver::rebuildOrderHeap()
{
    vec<Var> vs;
    for (Var v = 0; v < nVars(); v++)
        if (decision[v] && value(v) == l_Undef)
            vs.push(v);

    order_heap_CHB  .build(vs);
    order_heap_VSIDS.build(vs);
}


/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|  
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || propagate() != CRef_Undef)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    safeRemoveSatisfiedCompact(learnts_core, CORE);
    safeRemoveSatisfiedCompact(learnts_tier2, TIER2);
    safeRemoveSatisfiedCompact(learnts_local, LOCAL);

    if (remove_satisfied)        // Can be turned off.
        removeSatisfied(clauses);
    checkGarbage();
    rebuildOrderHeap();

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}


/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|  
|  Description:
|    Search for a model the specified number of conflicts. 
|  
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int& nof_conflicts)
{
    assert(ok);
    int         backtrack_level;
    int         lbd;
    vec<Lit>    learnt_clause;
    bool        cached = false;
    starts++;

    for (;;){
       if (decisionLevel() == 0) {
          importUnitClauses();
          if (!importClauses()) return l_False;
       }

       CRef confl = propagate();

       if (confl != CRef_Undef){
          // CONFLICT
          if (VSIDS){
             if (--timer == 0 && var_decay < 0.95) timer = 5000, var_decay += 0.01;
          }else
             if (step_size > min_step_size) step_size -= step_size_dec;

          conflicts++; nof_conflicts--;
          if (conflicts == 100000 && learnts_core.size() < 100) core_lbd_cut = 5;
          if (decisionLevel() == 0) return l_False;

          learnt_clause.clear();
          analyze(confl, learnt_clause, backtrack_level, lbd);

          if (exportClauseCallback != NULL) {
             exportClauseCallback(issuer, lbd, learnt_clause);
          }

          cancelUntil(backtrack_level);

          if (BMM && learnt_clause.size() <= 2)
          {
              bayesian_update(learnt_clause);
              for( int i=0; i<learnt_clause.size(); i++ ) {
                  Var v = var(learnt_clause[i]);
                  polarity[v] = (parameters[v].a > parameters[v].b) ? false : true;
                  activity_VSIDS[v] = activity_CHB[v] = max(parameters[v].a,parameters[v].b)/(parameters[v].a+parameters[v].b);
              }
          }
          lbd--;
          if (VSIDS){
             cached = false;
             conflicts_VSIDS++;
             lbd_queue.push(lbd);
             global_lbd_sum += (lbd > 50 ? 50 : lbd); }

          if (learnt_clause.size() == 1){
             uncheckedEnqueue(learnt_clause[0]);
          }else{
             CRef cr = ca.alloc(learnt_clause, true);
             ca[cr].set_lbd(lbd);
             if (lbd <= core_lbd_cut){
                learnts_core.push(cr);
                ca[cr].mark(CORE);
             }else if (lbd <= tier2_lbd_cut){
                learnts_tier2.push(cr);
                ca[cr].mark(TIER2);
                ca[cr].touched() = conflicts;
             }else{
                learnts_local.push(cr);
                claBumpActivity(ca[cr]); }
             attachClause(cr);
             uncheckedEnqueue(learnt_clause[0], cr);
          }
          if (drup_file){
#ifdef BIN_DRUP
             binDRUP('a', learnt_clause, drup_file);
#else
             for (int i = 0; i < learnt_clause.size(); i++)
                fprintf(drup_file, "%i ", (var(learnt_clause[i]) + 1) * (-2 * sign(learnt_clause[i]) + 1));
             fprintf(drup_file, "0\n");
#endif
          }

          if (VSIDS) varDecayActivity();
          claDecayActivity();

          /*if (--learntsize_adjust_cnt == 0){
            learntsize_adjust_confl *= learntsize_adjust_inc;
            learntsize_adjust_cnt    = (int)learntsize_adjust_confl;
            max_learnts             *= learntsize_inc;

            if (verbosity >= 1)
            printf("c | %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", 
            (int)conflicts, 
            (int)dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(), (int)clauses_literals, 
            (int)max_learnts, nLearnts(), (double)learnts_literals/nLearnts(), progressEstimate()*100);
            }*/

       }else{
          // NO CONFLICT
          bool restart = false;
          if (!VSIDS)
             restart = nof_conflicts <= 0;
          else if (!cached){
             restart = lbd_queue.full() && (lbd_queue.avg() * 0.8 > global_lbd_sum / conflicts_VSIDS);
             cached = true;
          }
          if (restart /*|| !withinBudget()*/){
             lbd_queue.clear();
             cached = false;
             // Reached bound on number of conflicts:
             progress_estimate = progressEstimate();
             cancelUntil(0);
             return l_Undef; }

          // Simplify the set of problem clauses:
          if (decisionLevel() == 0 && !simplify())
             return l_False;

          if (conflicts >= next_T2_reduce){
             next_T2_reduce = conflicts + 10000;
             reduceDB_Tier2(); }
          if (conflicts >= next_L_reduce){
             next_L_reduce = conflicts + 15000;
             reduceDB(); }

          Lit next = lit_Undef;
          while (decisionLevel() < assumptions.size()){
          // Perform user provided assumption:
          Lit p = assumptions[decisionLevel()];
          if (value(p) == l_True){
          // Dummy decision level:
          newDecisionLevel();
          }else if (value(p) == l_False){
          analyzeFinal(~p, conflict);
          return l_False;
          }else{
          next = p;
          break;
          }
          }

          if (next == lit_Undef){
             // New variable decision:
             decisions++;
             next = pickBranchLit();
             if (next == lit_Undef)
                // Model found:
                return l_True;
          }

          // Increase decision level and enqueue 'next'
          newDecisionLevel();
          uncheckedEnqueue(next);
       }
    }
}


double Solver::progressEstimate() const
{
   double  progress = 0;
   double  F = 1.0 / nVars();

   for (int i = 0; i <= decisionLevel(); i++){
      int beg = i == 0 ? 0 : trail_lim[i - 1];
      int end = i == decisionLevel() ? trail.size() : trail_lim[i];
      progress += pow(F, i) * (end - beg);
   }

   return progress / nVars();
}

/*
   Finite subsequences of the Luby-sequence:

0: 1
1: 1 1 2
2: 1 1 2 1 1 2 4
3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
...


 */

static double luby(double y, int x){

   // Find the finite subsequence that contains index 'x', and the
   // size of that subsequence:
   int size, seq;
   for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

   while (size-1 != x){
      size = (size-1)>>1;
      seq--;
      x = x % size;
   }

   return pow(y, seq);
}

template <typename T>
void Solver::bayesian_update(T& c)
{
    double coeff_product = 1.0;
    for( int j=0; j<c.size(); j++ )
    {
        Lit l = c[j];
        if ( value(l) == l_True ) return; // this clause is already satisfied
        if ( value(l) == l_False ) continue; // this literal is already falsified
        int v = var(l);
        int sgn = !sign(l);

        //if ( v >= nVars() ) printf("DBG: var %d is out of range %d\n", v, nVars());
        //if ( v >= parameters.size() || v >= updatedParams.size() ) printf("DBG: var %d (lit %d) is out of range (up.size: %d, p.size: %d)\n", v, toInt(l), updatedParams.size(), parameters.size()); 

        updatedParams[v].a = parameters[v].a + (!sgn);
        updatedParams[v].b = parameters[v].b + (sgn);

        double coeff = (sgn ? parameters[v].b : parameters[v].a) / (parameters[v].a + parameters[v].b);
        coeff_product *= coeff;
    }
    double normalization_constant = 1 - coeff_product;

    for( int j=0; j<c.size(); j++ )
    {
        Lit l = c[j];
        if ( value(l) == l_False ) continue; // this literal is already falsified
        int v = var(l);

        double sumP = parameters[v].a + parameters[v].b;
        double sumUP = updatedParams[v].a + updatedParams[v].b;

        double *p[2], *up[2];
        p[0] = &parameters[v].a;
        p[1] = &parameters[v].b;
        up[0] = &updatedParams[v].a;
        up[1] = &updatedParams[v].b;
        for( int k=0; k<=1; k++ )
        {
            double moment1 = *p[k] / sumP - coeff_product * *up[k] / sumUP;
            double moment2 = *p[k] * (*p[k] + 1) / ((sumP) * (sumP+1)) - coeff_product * *up[k] * (*up[k] + 1) / ((sumUP) * (sumUP+1));

            moment1 /= normalization_constant;
            moment2 /= normalization_constant;

            *p[k] = moment1 * (moment1 - moment2) / (moment2 - moment1*moment1);
        }
    }
}

void Solver::bayesian()
{
    BMM = 1;
    int epochs = 10;//bayesian_init_epochs;
    int n = nVars();

    updatedParams.growTo(n);

    for( int k=0; k<epochs; k++ )
    {
        for( int i=0; i<nClauses(); i++ )
        {
            Clause& c = ca[clauses[i]];
            bayesian_update(c);
        }

/*        for( int i=0; i<nLearnts(); i++ )
        {
            Clause& c = ca[learnts[i]];
            bayesian_update(c);
        }*/
    }

    for( int v=0; v<n; v++ )
    {
        polarity[v] = (parameters[v].a > parameters[v].b) ? false : true;
        activity_VSIDS[v] = activity_CHB[v] =
            max(parameters[v].a, parameters[v].b)/(parameters[v].a+parameters[v].b);
    }
}

void Solver::init_bayesian()
{

    int n = nVars();
    parameters.growTo(n);
    for( int i=0; i<n; i++ )
    {
        parameters[i].a = (double)rand() / RAND_MAX + 10;
        parameters[i].b = (double)rand() / RAND_MAX + 10;
    }

}

void Solver::jeroslow_wang_init(bool act, bool pol)
{
    if ( !act && !pol ) return;
    
    vec<double> cnt;
    int n = nVars();
    cnt.growTo(2 * n);
    for( int v=0; v<2*n; v++ )
        cnt[v] = 0.0;

    for( int i=0; i<nClauses(); i++ )
    {
        Clause& c = ca[clauses[i]];
        if ( c.size() > 50 ) continue;
        double sc = pow(2, -c.size());
        for( int j=0; j<c.size(); j++ )
        {
            Lit q = c[j];
            if ( sign(q) )
                cnt[var(q) + n] += sc;
            else
                cnt[var(q)] += sc;
        }
    }

    jw_pol = pol;
    jw_act = act;

    if ( jw_pol )
    {
        for( int v=0; v<n; v++ )
            polarity[v] = (cnt[v] > cnt[v+n]) ? false : true;
    }

    if ( jw_act )
    {
        int m = nClauses();
        if ( m == 0 ) m = 1;
        for( int v=0; v<n; v++ )
        {
            activity_CHB[v] = activity_VSIDS[v] = (cnt[v] + cnt[v+n]) / m;
        }
    }
}


//static bool switch_mode = false;
//static void SIGALRM_switch(int signum) { switch_mode = true; }

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_()
{
   //signal(SIGALRM, SIGALRM_switch);
   //alarm(2500);

   model.clear();
   conflict.clear();
   if (!ok) return l_False;

   solves++;

   max_learnts               = nClauses() * learntsize_factor;
   learntsize_adjust_confl   = learntsize_adjust_start_confl;
   learntsize_adjust_cnt     = (int)learntsize_adjust_confl;
   lbool   status            = l_Undef;

   if (verbosity >= 1){
      printf("c ============================[ Search Statistics ]==============================\n");
      printf("c | Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
      printf("c |           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
      printf("c ===============================================================================\n");
   }

   add_tmp.clear();

   //VSIDS = true;
   //int init = 10000;
   //while (status == l_Undef && init > 0 && withinBudget())
   //   status = search(init);
   //VSIDS = false;

   // Search:
   int curr_restarts = 0;
   while (status == l_Undef && withinBudget()){
      if (VSIDS){
         int weighted = INT32_MAX;
         status = search(weighted);
      }else{
         int nof_conflicts = luby(restart_inc, curr_restarts) * restart_first;
         curr_restarts++;
         status = search(nof_conflicts);
      }
//      if (!VSIDS && switch_mode){
//         VSIDS = true;
//         printf("c Switched to VSIDS.\n");
//         fflush(stdout);
//         picked.clear();
//         conflicted.clear();
//         almost_conflicted.clear();
//#ifdef ANTI_EXPLORATION
//         canceled.clear();
//#endif
//      }
   }

   if (verbosity >= 1)
      printf("c ===============================================================================\n");

#ifdef BIN_DRUP
   if (drup_file && status == l_False) binDRUP_flush(drup_file);
#endif

   if (status == l_True){
      // Extend & copy model:
      model.growTo(nVars());
      for (int i = 0; i < nVars(); i++) model[i] = value(i);
   }else if (status == l_False && conflict.size() == 0)
      ok = false;

   cancelUntil(0);
   return status;
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max)
{
   if (map.size() <= x || map[x] == -1){
      map.growTo(x+1, -1);
      map[x] = max++;
   }
   return map[x];
}


void Solver::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max)
{
   if (satisfied(c)) return;

   for (int i = 0; i < c.size(); i++)
      if (value(c[i]) != l_False)
         fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max)+1);
   fprintf(f, "0\n");
}


void Solver::toDimacs(const char *file, const vec<Lit>& assumps)
{
   FILE* f = fopen(file, "wr");
   if (f == NULL)
      fprintf(stderr, "could not open file %s\n", file), exit(1);
   toDimacs(f, assumps);
   fclose(f);
}


void Solver::toDimacs(FILE* f, const vec<Lit>& assumps)
{
   // Handle case when solver is in contradictory state:
   if (!ok){
      fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
      return; }

   vec<Var> map; Var max = 0;

   // Cannot use removeClauses here because it is not safe
   // to deallocate them at this point. Could be improved.
   int cnt = 0;
   for (int i = 0; i < clauses.size(); i++)
      if (!satisfied(ca[clauses[i]]))
         cnt++;

   for (int i = 0; i < clauses.size(); i++)
      if (!satisfied(ca[clauses[i]])){
         Clause& c = ca[clauses[i]];
         for (int j = 0; j < c.size(); j++)
            if (value(c[j]) != l_False)
               mapVar(var(c[j]), map, max);
      }

   // Assumptions are added as unit clauses:
   cnt += assumptions.size();

   fprintf(f, "p cnf %d %d\n", max, cnt);

   for (int i = 0; i < assumptions.size(); i++){
      assert(value(assumptions[i]) != l_False);
      fprintf(f, "%s%d 0\n", sign(assumptions[i]) ? "-" : "", mapVar(var(assumptions[i]), map, max)+1);
   }

   for (int i = 0; i < clauses.size(); i++)
      toDimacs(f, ca[clauses[i]], map, max);

   if (verbosity > 0)
      printf("c Wrote %d clauses with %d variables.\n", cnt, max);
}


//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to)
{
   // All watchers:
   //
   // for (int i = 0; i < watches.size(); i++)
   watches.cleanAll();
   watches_bin.cleanAll();
   for (int v = 0; v < nVars(); v++)
      for (int s = 0; s < 2; s++){
         Lit p = mkLit(v, s);
         // printf(" >>> RELOCING: %s%d\n", sign(p)?"-":"", var(p)+1);
         vec<Watcher>& ws = watches[p];
         for (int j = 0; j < ws.size(); j++)
            ca.reloc(ws[j].cref, to);
         vec<Watcher>& ws_bin = watches_bin[p];
         for (int j = 0; j < ws_bin.size(); j++)
            ca.reloc(ws_bin[j].cref, to);
      }

   // All reasons:
   //
   for (int i = 0; i < trail.size(); i++){
      Var v = var(trail[i]);

      if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
         ca.reloc(vardata[v].reason, to);
   }

   // All learnt:
   //
   for (int i = 0; i < learnts_core.size(); i++)
      ca.reloc(learnts_core[i], to);
   for (int i = 0; i < learnts_tier2.size(); i++)
      ca.reloc(learnts_tier2[i], to);
   for (int i = 0; i < learnts_local.size(); i++)
      ca.reloc(learnts_local[i], to);

   // All original:
   //
   int i, j;
   for (i = j = 0; i < clauses.size(); i++)
      if (ca[clauses[i]].mark() != 1){
         ca.reloc(clauses[i], to);
         clauses[j++] = clauses[i]; }
   clauses.shrink(i - j);
}


void Solver::garbageCollect()
{
   // Initialize the next region to a size corresponding to the estimated utilization degree. This
   // is not precise but should avoid some unnecessary reallocations for the new region:
   ClauseAllocator to(ca.size() - ca.wasted()); 

   relocAll(to);
   if (verbosity >= 2)
      printf("c |  Garbage collection:   %12d bytes => %12d bytes             |\n", 
            ca.size()*ClauseAllocator::Unit_Size, to.size()*ClauseAllocator::Unit_Size);
   to.moveTo(ca);
}
