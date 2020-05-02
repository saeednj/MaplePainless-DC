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

#include <atomic>
#include <stdlib.h>

#include "../clauses/ClauseExchange.h"

using namespace std;

/// Class in charge of the management of shared clauses.
class ClauseManager
{
public:
   /// Init the clause manager.
   static void initClauseManager()
   {
   }

   /// Alloc a new shared clause.
   static ClauseExchange * allocClause(int size)
   {
      ClauseExchange * ptr;
      ptr = (ClauseExchange *) malloc(sizeof(ClauseExchange) + sizeof(int) * size);

      ptr->size   = size;
      ptr->nbRefs = 1;

      return ptr;
   }

   /// Increase the number of references to a shared clause.
   static void increaseClause(ClauseExchange * cls, int refs = 1)
   {
     /*if(cls==NULL){
       printf("increaseClause:cls NULL.\n");
       abort();
       }*/

      cls->nbRefs += refs; // atomic adition
   }

   /// Release a shared clauses, delete it if needed.
   static void releaseClause(ClauseExchange * cls)
   {
     /*if(cls==NULL){
       printf("releaseClause:cls NULL.\n");
       abort();
       }*/

      int oldValue = cls->nbRefs.fetch_sub(1); // atomic decrementation

      if (oldValue - 1 <= 0) {
         // Only the last thread should execute this code
	/*if(oldValue-1<=-1){
         fprintf(stderr,"\x1B[31mSomeone tried to free clause %p more than once.\x1B[0m\n",cls);
	 }*/
         free(cls);
      }
   }
   
   /// Join the clause manager.
   static void joinClauseManager()
   {
   }
};
