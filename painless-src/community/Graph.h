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

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

class Graph {
public:
   unsigned      nNodes;         ///< Number of nodes.
   unsigned long nEdges;         ///< Number of edges.
   double        totalWeight;    ///< Total weight.

   vector<unsigned long> nodes;    ///< Vector of nodes.
   vector<unsigned>      edges;    ///< Vector of edges.
   vector<double>        weights;  ///< Vector of edge weight.
   
   // CONSTRUCTORS
   
   Graph();
   Graph(vector< map<int, double> > & adjacencyList);
   Graph(int nNodes, int nEdges, double totalWeight,
         int * nodes, int * edges, double * weights);

   // binary file format is
   // 4 bytes for the number of nodes in the graph
   // 8*(nb_nodes) bytes for the cumulative degree for each node:
   //    deg(0)=nodes[0]
   //    deg(k)=nodes[k]-nodes[k-1]
   // 4*(sum_nodes) bytes for the edges
   // IF WEIGHTED 4*(sum_nodes) bytes for the weights in a separate file
   Graph(char *filename, char *filename_w, int type);

   // TO DISPLAY GRAPH 
   void displayAdjacencyList(FILE * output) const;
   void displayEdgeList(FILE * output) const;
   void displayReverseEdgeList(FILE * output) const;
   void displayBinary(const char * outfile) const;
   bool checkSymmetry(FILE * output) const;

   // return the number of neighbors (degree) of the node
   inline unsigned nNeighbors(unsigned node) const;

   // return the number of self loops of the node
   inline double nSelfloops(unsigned node) const;

   // return the weighted degree of the node
   inline double weightedDegree(unsigned node) const;

   // return pointers to the first neighbor and first weight of the node
   inline pair<vector<unsigned>::const_iterator, vector<double>::const_iterator >
      neighbors(unsigned node) const;
};

inline unsigned Graph::nNeighbors(unsigned node) const
{
   assert(node < nNodes);

   if (node == 0)
      return nodes[0];

   return nodes[node] - nodes[node-1];
}

inline double Graph::nSelfloops(unsigned node) const
{
   assert(node < nNodes);

   auto p = neighbors(node);
   for (unsigned i = 0; i < nNeighbors(node); i++) {
      if (*(p.first+i) == node) {
         if (weights.size() != 0)
            return (double) *(p.second + i);
         return 1.0;
      }
   }

   return 0.0;
}

inline double Graph::weightedDegree(unsigned node) const
{
   assert(node < nNodes);

   if (weights.size() == 0)
      return (double)nNeighbors(node);

   auto p = neighbors(node);
   double res = 0;
   for (unsigned i = 0; i < nNeighbors(node); i++) {
      res += (double) *(p.second + i);
   }

   return res;
}

inline pair<vector<unsigned>::const_iterator, vector<double>::const_iterator >
Graph::neighbors(unsigned node) const
{
   assert(node < nNodes);

   if (node == 0)
      return make_pair(edges.begin(), weights.begin());

   if (weights.size() != 0)
      return make_pair(edges.begin() + nodes[node-1], weights.begin() +
                       nodes[node-1]);

   return make_pair(edges.begin() + nodes[node-1], weights.begin());
}
