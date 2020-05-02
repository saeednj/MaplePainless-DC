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


#include "../working/SimpleWorkingQueue.h"


SimpleWorkingQueue::SimpleWorkingQueue()
{
}

SimpleWorkingQueue::~SimpleWorkingQueue()
{
}

void
SimpleWorkingQueue::pushJob(vector<int> & guidingPath)
{
   queueLock.lock();
   queue.push_back(guidingPath);
   queueLock.unlock();
}
   
bool
SimpleWorkingQueue::popJob(vector<int> & outGuidingPath)
{
   queueLock.lock();

   if (queue.empty()) {
      queueLock.unlock();
      return false;
   }

   outGuidingPath = queue.front();
   queue.pop_front();

   queueLock.unlock();

   return true;
}
   
bool
SimpleWorkingQueue::empty()
{
   queueLock.lock();
   bool ret = queue.empty();
   queueLock.unlock();

   return ret;
}
