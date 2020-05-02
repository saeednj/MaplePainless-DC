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

#include "../clauses/ClauseExchange.h"

#include <bitset>
#include <vector>

using namespace std;

//#define NUM_BITS 268435399 // 32MB
#define NUM_BITS 26843543 // 3,2MB

/// Bloom filter for clauses
class ClauseFilter
{
public:
   /// Constructor.
	ClauseFilter();

   /// Destructor.
	virtual ~ClauseFilter();
	 
   /// Try to add the clause to the filter.
	bool registerClause(ClauseExchange * cls);
	
   /// Clear the filter, i.e., return to its initial state.
	void clear();

protected:
   /// Used to store the hashes of the shared clauses.
	bitset<NUM_BITS>* s1;

   /// The hash function.
	size_t commutativeHashFunction(ClauseExchange * cls, int which);
};
