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

#include <limits>

#include "../community/Community.h"
#include "../utils/Logger.h"

using namespace std;

Community::Community(char * filename, char * filename_w, int type,
                     double precision)
{
   g    = Graph(filename, filename_w, type);
   size = g.nNodes;

   neigh_weight.resize(size,-1);
   neigh_pos.resize(size);
   neigh_last = 0;

   n2c.resize(size);
   in.resize(size);
   tot.resize(size);

   for (int i = 0; i < size; i++) {
      n2c[i] = i;
      tot[i] = g.weightedDegree(i);
      in[i]  = g.nSelfloops(i);
   }

   precision = precision;
}

Community::Community(Graph gc, double precision)
{
   g    = gc;
   size = g.nNodes;

   neigh_weight.resize(size,-1);
   neigh_pos.resize(size);
   neigh_last = 0;

   n2c.resize(size);
   in.resize(size);
   tot.resize(size);

   for (int i = 0; i < size; i++) {
      n2c[i] = i;
      in[i]  = g.nSelfloops(i);
      tot[i] = g.weightedDegree(i);
   }

   this->precision = precision;
}

void
Community::init_partition(char * filename)
{
   ifstream finput;
   finput.open(filename,fstream::in);

   // read partition
   while (!finput.eof()) {
      unsigned node, comm;
      finput >> node >> comm;

      if (finput) {
         int old_comm = n2c[node];
         neigh_comm(node);

         remove(node, old_comm, neigh_weight[old_comm]);

         unsigned i = 0;
         for (i = 0 ; i < neigh_last ; i++) {
            unsigned best_comm = neigh_pos[i];
            double best_nblinks     = neigh_weight[neigh_pos[i]];
            if (best_comm == comm) {
               insert(node, best_comm, best_nblinks);
               break;
            }
         }
         if (i == neigh_last)
            insert(node, comm, 0);
      }
   }
   finput.close();
}

void
Community::display()
{
   for (int i = 0; i < size; i++)
      cerr << " " << i << "/" << n2c[i] << "/" << in[i] << "/" << tot[i] ;
   cerr << endl;
}

double
Community::modularity()
{
   double q  = 0.;
   double m2 = (double)g.totalWeight;

   for (int i = 0; i < size; i++) {
      if (tot[i] > 0)
         q += (double)in[i] / m2 - ((double)tot[i] / m2) * ((double)tot[i] / m2);
   }

   return q;
}

void
Community::neigh_comm(unsigned node)
{
   for (unsigned i = 0; i < neigh_last; i++)
      neigh_weight[neigh_pos[i]] = -1;
   neigh_last = 0;

   auto p = g.neighbors(node);

   unsigned deg = g.nNeighbors(node);

   neigh_pos[0]               = n2c[node];
   neigh_weight[neigh_pos[0]] = 0;
   neigh_last                 = 1;

   for (unsigned i = 0; i < deg; i++) {
      unsigned neigh       = *(p.first+i);
      unsigned _neigh_comm = n2c[neigh];
      double neigh_w       = (g.weights.size() == 0) ? 1. : *(p.second + i);

      if (neigh != node) {
         if (neigh_weight[_neigh_comm] == -1) {
            neigh_weight[_neigh_comm] = 0.;
            neigh_pos[neigh_last++]   = _neigh_comm;
         }
         neigh_weight[_neigh_comm] += neigh_w;
      }
   }
}

void
Community::partition2graph()
{
   vector<int> renumber(size, -1);
   for (int node = 0; node < size; node++) {
      renumber[n2c[node]]++;
   }

   int final = 0;
   for (int i = 0; i < size; i++)
      if (renumber[i] != -1)
         renumber[i]=final++;

   for (int i = 0; i < size; i++) {
      auto p = g.neighbors(i);

      int deg = g.nNeighbors(i);
      for (int j = 0; j < deg; j++) {
         int neigh = *(p.first+j);
         cout << renumber[n2c[i]] << " " << renumber[n2c[neigh]] << endl;
      }
   }
}

void
Community::display_partition() {
   vector<int> renumber(size, -1);
   for (int node = 0; node < size; node++) {
      renumber[n2c[node]]++;
   }

   int final = 0;
   for (int i = 0; i < size; i++)
      if (renumber[i] != -1)
         renumber[i] = final++;

   for (int i = 0; i < size; i++)
      cout << i << " " << renumber[n2c[i]] << endl;
}

vector< pair<int, int> >
Community::get_partition()
{
   vector< pair<int, int> > ret(0);
   vector<int> renumber(size, -1);
   for (int node = 0; node < size; node++) {
      renumber[n2c[node]]++;
   }

   int final = 0;
   for (int i = 0; i < size; i++)
      if (renumber[i] != -1)
         renumber[i] = final++;

   for (int i = 0; i < size; i++)
      ret.push_back(make_pair(i, renumber[n2c[i]]));
   //cout << i << " " << renumber[n2c[i]] << endl;

   return ret;
}

Graph
Community::partition2graph_binary() {
   // Renumber communities
   vector<int> renumber(size, -1);
   for (int node = 0; node < size; node++) {
      renumber[n2c[node]]++;
   }

   int final = 0;
   for (int i = 0; i < size; i++)
      if (renumber[i] != -1)
         renumber[i] = final++;

   // Compute communities
   vector<vector<int> > comm_nodes(final);
   for (int node = 0; node < size; node++) {
      comm_nodes[renumber[n2c[node]]].push_back(node);
   }

   // Compute weighted graph
   Graph g2;
   g2.nNodes = comm_nodes.size();
   g2.nodes.resize(comm_nodes.size());

   int comm_deg = comm_nodes.size();
   for (int comm = 0; comm < comm_deg; comm++) {
      map<int,double> m;
      map<int,double>::iterator it;

      int comm_size = comm_nodes[comm].size();
      for (int node = 0; node < comm_size; node++) {
         auto p = g.neighbors(comm_nodes[comm][node]);
         int deg = g.nNeighbors(comm_nodes[comm][node]);

         for (int i = 0; i < deg; i++) {
            int neigh           = *(p.first+i);
            int neigh_comm      = renumber[n2c[neigh]];
            double neigh_weight = (g.weights.size() == 0) ? 1. : *(p.second+i);

            it = m.find(neigh_comm);
            if (it == m.end())
               m.insert(make_pair(neigh_comm, neigh_weight));
            else
               it->second += neigh_weight;
         }
      }
      g2.nodes[comm] = (comm == 0) ? m.size() : g2.nodes[comm - 1] + m.size();
      g2.nEdges += m.size();

      for (it = m.begin(); it != m.end(); it++) {
         g2.totalWeight  += it->second;
         g2.edges.push_back(it->first);
         g2.weights.push_back(it->second);
      }
   }

   return g2;
}

bool
Community::one_level()
{
   int nb_moves;
   bool improvement = false ;
   int nb_pass_done = 0;
   double new_mod   = modularity();
   double cur_mod   = new_mod;

   vector<int> random_order(size);
   for (int i = 0; i < size; i++)
      random_order[i] = i;
   for (int i = 0; i < size - 1; i++) {
      int rand_pos           = rand() % (size - i) + i;
      int tmp                = random_order[i];
      random_order[i]        = random_order[rand_pos];
      random_order[rand_pos] = tmp;
   }

   // repeat while
   //   there is an improvement of modularity
   //   or there is an improvement of modularity greater than a given epsilon
   //   or a predefined number of pass have been done
   do {
      cur_mod  = new_mod;
      nb_moves = 0;
      nb_pass_done++;

      // for each node: remove the node from its community and insert it in the
      // best community
      for (int node_tmp = 0; node_tmp < size; node_tmp++) {
         //      int node = node_tmp;
         int node        = random_order[node_tmp];
         int node_comm   = n2c[node];
         double w_degree = g.weightedDegree(node);

         // computation of all neighboring communities of current node
         neigh_comm(node);
         // remove node from its current community
         remove(node, node_comm, neigh_weight[node_comm]);

         // compute the nearest community for node
         // default choice for future insertion is the former community
         int best_comm        = node_comm;
         double best_nblinks  = 0.;
         double best_increase = 0.;
         for (unsigned i = 0; i < neigh_last; i++) {
            double increase = modularity_gain(node, neigh_pos[i],
                  neigh_weight[neigh_pos[i]], w_degree);
            if (increase > best_increase) {
               best_comm     = neigh_pos[i];
               best_nblinks  = neigh_weight[neigh_pos[i]];
               best_increase = increase;
            }
         }

         // insert node in the nearest community
         insert(node, best_comm, best_nblinks);

         if (best_comm != node_comm)
            nb_moves++;
      }

      double total_tot = 0;
      double total_in  = 0;
      for (unsigned i = 0; i < tot.size(); i++) {
         total_tot += tot[i];
         total_in  += in[i];
      }

      new_mod = modularity();
      if (nb_moves > 0)
         improvement = true;

      log(0, "nb_pass %d mod=%f\n", nb_pass_done, new_mod);
   } while (nb_moves > 0 && new_mod - cur_mod > precision);

            return improvement;
}
