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

#include <stdlib.h>
#include <stdio.h>

#include "../clauses/ClauseFilter.h"

using namespace std;

#define NUM_PRIMES 12

static unsigned const int primes [] = { 2038072819, 2038073287, 2038073761,
                                        2038074317, 2038072823, 2038073321,
                                        2038073767, 2038074319, 2038072847,
                                        2038073341, 2038073789, 2038074329 };

size_t
ClauseFilter::commutativeHashFunction(ClauseExchange * cls, int which)
{
	size_t res = 0;

	for (size_t j = 0; j < cls->size; j++) {
		int lit = cls->lits[j];

		res ^= lit * primes[abs((which * lit) % NUM_PRIMES)];
	}

	return res % NUM_BITS;
}

ClauseFilter::ClauseFilter()
{
	s1 = new bitset<NUM_BITS>();
}

ClauseFilter::~ClauseFilter()
{
	delete s1;
}

bool
ClauseFilter::registerClause(ClauseExchange * cls)
{
	// unit clauses always get in
	if (cls->size == 1)
		return true;

	size_t h1 = commutativeHashFunction(cls, 1);
	size_t h2 = commutativeHashFunction(cls, 2);
	size_t h3 = commutativeHashFunction(cls, 3);
	size_t h4 = commutativeHashFunction(cls, 4);

	if (s1->test(h1) && s1->test(h2) && s1->test(h3) && s1->test(h4))
		return false;

	s1->set(h1, true);
	s1->set(h2, true);
	s1->set(h3, true);
	s1->set(h4, true);
	
   return true;
}

void
ClauseFilter::clear()
{
	s1->reset();
}
