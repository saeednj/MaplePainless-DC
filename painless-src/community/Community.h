//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large
// networks" Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E.
// Lefebvre
//
// This file is part of Louvain algorithm.
// 
// Louvain algorithm is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
// 
// Louvain algorithm is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
// 
// You should have received a copy of the GNU General Public License along with
// Louvain algorithm. If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time     : February 2008
//-----------------------------------------------------------------------------

#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <sys/mman.h>
#include <assert.h>

#include "Graph.h"

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

class Community {
public:
   vector<double> neigh_weight;
   vector<unsigned> neigh_pos;
   unsigned neigh_last;

   Graph g; // network to compute communities for
   int size; // nummber of nodes in the network and size of all vectors
   vector<int> n2c; // community to which each node belongs
   vector<double> in,tot; // used to compute the modularity participation of
                          // each community

   // number of pass for one level computation
   // if -1, compute as many pass as needed to increase modularity
   int nb_pass;

   // a new pass is computed if the last one has generated an increase 
   // greater than precision
   // if 0. even a minor increase is enough to go for one more pass
   double precision;

   // constructors:
   // reads graph from file using graph constructor
   // type defined the weighted/unweighted status of the graph file
   Community (char *filename, char *filename_w, int type, double precision);
   // copy graph
   Community (Graph g, double precision);

   // initiliazes the partition with something else than all nodes alone
   void init_partition(char *filename_part);

   // display the community of each node
   void display();

   // remove the node from its current community with which it has dnodecomm
   // links
   inline void remove(int node, int comm, double dnodecomm);

   // insert the node in comm with which it shares dnodecomm links
   inline void insert(int node, int comm, double dnodecomm);

   // compute the gain of modularity if node where inserted in comm
   // given that node has dnodecomm links to comm.  The formula is:
   // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
   // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
   // where In(comm)    = number of half-links strictly inside comm
   //       Tot(comm)   = number of half-links inside or outside comm
   //                     (sum(degrees))
   //       d(node,com) = number of links from node to comm
   //       deg(node)   = node degree
   //       m           = number of links
   inline double modularity_gain(int node, int comm, double dnodecomm,
                                 double w_degree);

   // compute the set of neighboring communities of node
   // for each community, gives the number of links from node to comm
   void neigh_comm(unsigned node);

   // compute the modularity of the current partition
   double modularity();

   // displays the graph of communities as computed by one_level
   void partition2graph();
   // displays the current partition (with communities renumbered from 0 to k-1)
   void display_partition();

   vector< pair<int, int> > get_partition();

   // generates the binary graph of communities as computed by one_level
   Graph partition2graph_binary();

   // compute communities of the graph for one level
   // return true if some nodes have been moved
   bool one_level();
};

inline void Community::remove(int node, int comm, double dnodecomm)
{
   assert(node >= 0 && node < size);

   tot[comm] -= g.weightedDegree(node);
   in[comm]  -= 2*dnodecomm + g.nSelfloops(node);
   n2c[node]  = -1;
}

inline void Community::insert(int node, int comm, double dnodecomm)
{
   assert(node >= 0 && node < size);

   tot[comm] += g.weightedDegree(node);
   in[comm]  += 2*dnodecomm + g.nSelfloops(node);
   n2c[node]  = comm;
}

inline double Community::modularity_gain(int node, int comm, double dnodecomm,
                                         double w_degree)
{
   assert(node >= 0 && node < size);

   double totc = (double)tot[comm];
   double degc = (double)w_degree;
   double m2   = (double)g.totalWeight;
   double dnc  = (double)dnodecomm;

   return (dnc - totc*degc/m2);
}
