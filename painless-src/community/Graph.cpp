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

#include "Graph.h"

using namespace std;

Graph::Graph()
{
   nNodes      = 0;
   nEdges      = 0;
   totalWeight = 0;
}

Graph::Graph(vector< map<int, double> > & adjacencyList)
{
   nNodes = adjacencyList.size();
   nodes.resize(nNodes);

   long long tot = 0;
   int count     = 0;

   for(unsigned i = 0; i < adjacencyList.size(); i++){
      tot += adjacencyList[i].size();
      nodes[i] = tot;
   }

   nEdges = tot;
   edges.resize(nEdges);

   weights.resize(nEdges);
   totalWeight = 0;
   count = 0;
   
  map<int, double>::iterator it;
   for(size_t i = 0; i < adjacencyList.size(); i++){
      for(it = adjacencyList[i].begin(); it != adjacencyList[i].end(); it++){
         edges[count]   = it->first;
         weights[count] = it->second;
         totalWeight   += it->second;
         count++;
      }
   }
}

Graph::Graph(char * filename, char * filename_w, int type)
{
   ifstream finput;
   finput.open(filename, fstream::in | fstream::binary);

   // Read number of nodes on 4 bytes
   finput.read((char *)&nNodes, 4);
   assert(finput.rdstate() == ios::goodbit);

   // Read cumulative degree sequence: 8 bytes for each node
   // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
   nodes.resize(nNodes);
   finput.read((char *)&nodes[0], nNodes*8);

   // Read edges: 4 bytes for each edge (each edge is counted twice)
   nEdges = nodes[nNodes-1];
   edges.resize(nEdges);
   finput.read((char *)(&edges[0]), (long)nEdges*8);  

   // IF WEIGHTED : read weights: 4 bytes for each edge (each edge is counted
   // twice)
   weights.resize(0);
   totalWeight=0;
   if (type == WEIGHTED) {
      ifstream finput_w;
      finput_w.open(filename_w, fstream::in | fstream::binary);
      weights.resize(nEdges);
      finput_w.read((char *)&weights[0], (long)nEdges*4);  
   }

   // Compute total weight
   for (unsigned i = 0; i < nNodes; i++) {
      totalWeight += (double)weightedDegree(i);
   }
}


void
Graph::displayAdjacencyList(FILE * output) const
{
   for (unsigned node = 0; node < nNodes; node++) {
      auto p = neighbors(node);
      fprintf(stdout, "%u:", node);
      
      for (unsigned i = 0; i < nNeighbors(node); i++) {
         if (weights.size() != 0) {
            fprintf(stdout, " (%u %f)", *(p.first+i), *(p.second+i));
         } else {
            fprintf(stdout, " %u", *(p.first+i));
         }
      }
      fprintf(output, "\n");
   }
}

void
Graph::displayEdgeList(FILE * output) const
{
   fprintf(output, "%u %lu\n", nNodes, nEdges);
   for (unsigned node = 0; node < nNodes; node++) {
      auto p = neighbors(node);
      for (unsigned i = 0; i < nNeighbors(node); i++) {
         fprintf(output, "%u %u\n", node, *(p.first+i));
      }
   }
}

void
Graph::displayReverseEdgeList(FILE * output) const
{
   for (unsigned node = 0; node < nNodes; node++) {
      auto p = neighbors(node);
      for (unsigned i = 0 ; i < nNeighbors(node); i++) {
         if (node > *(p.first+i)) {
            if (weights.size() != 0) {
               // cout << *(p.first+i) << " " << node << " " << *(p.second+i);
               // cout << endl;
               fprintf(output, "%u %u %f\n", *(p.first+i), node, *(p.second+i));
            }
            else
               // cout << *(p.first+i) << " " << node << endl;
               fprintf(output, "%u %u\n", *(p.first+i), node);
         }
      }
   }
}

bool
Graph::checkSymmetry(FILE * output) const
{
   int error=0;
   for (unsigned node = 0; node < nNodes; node++) {
      auto p = neighbors(node);
      for (unsigned i = 0; i<nNeighbors(node); i++) {
         unsigned neigh = *(p.first+i);
         double weight  = *(p.second+i);

         auto p_neigh = neighbors(neigh);
         for (unsigned j=0 ; j<nNeighbors(neigh) ; j++) {
            unsigned neigh_neigh = *(p_neigh.first+j);
            double neigh_weight       = *(p_neigh.second+j);

            if (node == neigh_neigh && weight != neigh_weight) {
               fprintf(output, "%u %u %f %f\n", node, neigh, weight,
                       neigh_weight);
               if (error++ == 10)
                  exit(0);
            }
         }
      }
   }
   return (error == 0);
}

void
Graph::displayBinary(const char * outfile) const
{
   ofstream foutput;
   foutput.open(outfile ,fstream::out | fstream::binary);

   foutput.write((char *)(&nNodes), 4);
   foutput.write((char *)(&nodes[0]), 4*nNodes);
   foutput.write((char *)(&edges[0]), 8*nEdges);
}
