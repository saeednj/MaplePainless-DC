#include "../community/Community.h"
#include "../community/Graph.h"
#include "../community/GraphAlgorithms.h"
#include "../utils/Logger.h"

vector<map<int, double> >
renumber(const vector< map<int, double> > & graph, vector<int> & r2o,
         vector<int> & o2r)
{
   r2o.clear();
   o2r.clear();
   vector< map<int, double> > ret;
   o2r.resize(graph.size());
   fill(o2r.begin(), o2r.end(), -1);
   for(size_t i = 0; i < graph.size(); i++) {
      if(graph[i].size() > 0){
         o2r[i] = r2o.size();
         r2o.push_back((int)i);
         //ret.push_back(graph[i]);
      }
   }

   int ns = r2o.size();
   ret.resize(ns);
   for(int i = 0; i < ns; i++){
      int original = r2o[i];
      for(map<int, double>::const_iterator it = graph[original].begin();
            it != graph[original].end(); it++)
         ret[i][o2r[it->first]] = it->second;
   }
   return ret;
}

void louvain(const Graph & graph, double precision, const vector<int> & r2o,
             vector<vector<int> > & area2nodes, vector<int> & nodes2area)
{
   srand(20150317);

   //Community c(Graph(graph), precision);
   Community c(graph, precision);

   Graph g;
   bool improvement = true;

   vector<pair<int, int> > partitions(0);

   do {
      improvement = c.one_level(); // segmentation fault

      // partition sequences
      vector< pair<int, int> > prt = c.get_partition();
      for(auto it = prt.begin(); it != prt.end(); it++) {
         partitions.push_back(*it);
      }

      g = c.partition2graph_binary();
      c = Community(g, precision);
   } while(improvement);

   //calc final partition
   vector< vector<int> > levels;
   int l = -1;
   for(int i = 0; i < partitions.size(); i++) {
      if(partitions[i].first == 0) {
         l++;
         levels.resize(l+1);
      }
      levels[l].push_back(partitions[i].second);
   }

   nodes2area.resize(levels[0].size());
   for(int i = 0; i < levels[0].size(); i++)
      nodes2area[i] = i;

   int mx_level = (int)levels.size();
   for(l = 0; l < mx_level; l++)
      for(int node = 0; node < levels[0].size(); node++)
         nodes2area[node] = levels[l][nodes2area[node]];

   area2nodes.resize(levels[mx_level-1].size());
   for (int node = 0; node < levels[0].size(); node++) {
      area2nodes[nodes2area[node]].push_back(r2o[node]);
      //area2nodes[nodes2area[node]].push_back(node);
   }
}

static int
countbit(int n)
{
   return (n == 0) ? 0 : (1 + countbit(n & (n - 1)));
}

void
pageRank(const Graph & graph, double * page_rank)
{
   int cnt      = 0;
   double ratio = 0.85;
   double *s    = new double[graph.nNodes];

   for (unsigned i = 0; i < graph.nNodes; i++) {
      if (graph.nNeighbors(i) > 0) {
         cnt++;
         page_rank[i] = 1;
      } else {
         page_rank[i] = 0;
      }
   }

   for (int step = 0; step < 100; step++) {
      for (unsigned i = 0; i < graph.nNodes; i++) {
         if (graph.nNeighbors(i) > 0) {
            s[i] = (double)(1.0 - ratio) / cnt;
         }

         auto p = graph.neighbors(i);
         for (unsigned j = 0 ; j < graph.nNeighbors(i) ; j++) {
            s[i] += ratio * page_rank[*(p.first+j)] /
               graph.nNeighbors(*(p.first+j));
         }
      }
      
      for (unsigned i = 0; i < graph.nNodes; i++) {
         if (graph.nNeighbors(i) > 0) {
            page_rank[i] = s[i];
         }
      }

      double maxp = 0;

      for (unsigned i = 0; i < graph.nNodes; i++) {
         if (page_rank[i] > maxp) {
            maxp = page_rank[i];
         }
      }

      if (maxp < 1e-10) {
         break;
      }

      for (unsigned i = 0; i < graph.nNodes; i++) {
         page_rank[i] /= maxp;
      }
   }

   delete[] s;
}

void
weightedPageRank(const Graph & graph, double * page_rank)
{
   int cnt      = 0;
   double ratio = 0.85;
   double *s    = new double[graph.nNodes];

   for (unsigned i = 0; i < graph.nNodes; i++) {
      if (graph.nNeighbors(i) > 0) {
         cnt++;
         page_rank[i] = 1;
      } else {
         page_rank[i] = 0;
      }
   }

   for (int step = 0; step < 50; step++) {
      for (unsigned i = 0; i < graph.nNodes; i++) {
         if (graph.nNeighbors(i) > 0) {
            s[i] = (double)(1.0 - ratio) / cnt;
         }

         auto p = graph.neighbors(i);
         for (unsigned j = 0 ; j < graph.nNeighbors(i) ; j++) {
            s[i] += ratio * page_rank[*(p.first+j)] * *(p.second+j) /
               graph.weightedDegree(*(p.first+j));
         }
      }

      for (unsigned i = 0; i < graph.nNodes; i++) {
         if (graph.nNeighbors(i) > 0) {
            page_rank[i] = s[i];
         }
      }

      double maxp = 0;

      for (unsigned i = 0; i < graph.nNodes; i++) {
         if (page_rank[i] > maxp) {
            maxp = page_rank[i];
         }
      }

      if (maxp < 1e-10) {
         break;
      }

      for (unsigned i = 0; i < graph.nNodes; i++) {
         page_rank[i] /= maxp;
      }
   }

   delete[] s;
}

static vector<int>
structure_hole_min_max_faster(int c, int size, double precision, double alpha,
                              double * beta, double * page_rank,
                              vector<int> target_communities, vector<int> area,
                              const Graph & graph)
{
   int heapsize          = 0;
   double ** influential = new double * [graph.nNodes];
   double ** sh          = new double * [graph.nNodes];
   int * pos             = new int      [graph.nNodes];
   double * tp           = new double   [graph.nNodes];
   int *heap             = new int      [graph.nNodes];

   for (unsigned i = 0; i < graph.nNodes; i++) {
      pos[i]         = -1;
      tp[i]          = 0;

      influential[i] = new double[c];
      for (int j = 0; j < c; j++) {
         influential[i][j] = 0.0;
      }

      sh[i] = new double[two(c)];
      for (int j = 0; j < two(c); j++) {
         sh[i][j] = 0.0;
      }
   }

   for (int k = 0; k < c; k++) {
      for (unsigned i = 0; i < graph.nNodes; i++) {
         if ((area[i] & target_communities[k]) == target_communities[k]) {
            influential[i][k] = page_rank[i];

            if (influential[i][k] > tp[i] + precision) {
               tp[i]            = influential[i][k];
               pos[i]           = heapsize-1;
               heap[heapsize++] = i;
               for (int i = heapsize-1, j = (i-1) >> 1; i > 0;
                    i = j, j = (i-1) >> 1) {

                  if (tp[heap[i]] > tp[heap[j]]) {
                     swap(heap[i],heap[j]);
                  } else {
                     break;
                  }
               }
            }
         }
      }
   }

   double *exp_inf = new double[c];
   while (heapsize > 0) {
      int idx = heap[0];
      tp[idx]  = 0;
      pos[idx] = -1;
      heap[0]  = heap[--heapsize];
      if (heapsize > 0) {
         int i = 0, key = heap[0];
         double tmp = tp[key];
         for (int j = (i << 1) + 1; j < heapsize; i = j, j = (i << 1) + 1) {
            if (j + 1 < heapsize && tp[heap[j+1]] > tp[heap[j]]) {
               j++;
            }
            if (tp[heap[j]] <= tmp)
               break;

            heap[i]      = heap[j];
            pos[heap[i]] = i;
         }
         heap[i]  = key;
         pos[key] = i;
      }

      double * gset = sh[idx], * gs = influential[idx];
      gset[0] = 1e100;
      for (int k = 0, set = 1; set < two(c); set++) {
         if (!contain(set,k)) {
            k++;
         }
         gset[set] = min(gs[k], gset[set ^ two(k)]);
      }
      for (int i = 0; i < c; i++) {
         exp_inf[i] = 0;
      }
      for (int set = 0; set < two(c); set++) {
         double tmp = beta[set] * gset[set];
         if (tmp < precision)
            continue;
         for (int i = 0; i < c; i++) {
            if (contain(set, i)) {
               checkmax(exp_inf[i], tmp);
            }
         }
      }

      auto p = graph.neighbors(idx);
      for (unsigned i = 0; i < graph.nNeighbors(idx); i++) {
         unsigned other = *(p.first+i);
         for (int i = 0; i < c; i++) {
            double new_inf = alpha * gs[i] + exp_inf[i];
            if (new_inf > influential[other][i] + precision) {
               influential[other][i] = new_inf;
               if (new_inf <= tp[other] + precision)
                  continue;
               tp[other] = new_inf;
               if (pos[other] < 0) {
                  pos[other]       = heapsize;
                  heap[heapsize++] = other;
               }
               int i = pos[other];
               for (int j = (i-1) >> 1; i > 0; i = j, j = (i-1) >> 1) {
                  if (tp[heap[j]] >= new_inf)
                     break;
                  heap[i]      = heap[j];
                  pos[heap[i]] = i;
               }
               heap[i]    = other;
               pos[other] = i;
            }
         }
      }
   }
   delete[] exp_inf;

   int qsize = 0;
   pair<double, int> * q = new pair<double, int> [graph.nNodes];
   for (unsigned i = 0; i < graph.nNodes; i++)
   {
      double s2 = 0;
      for (int k = 0; k < two(c); k++) {
         checkmax(s2, beta[k] * sh[i][k]);
      }
      double s3 = 0;
      for (int j = 0; j < c; j++) {
         s3 += influential[i][j];
      }
      double weight = (int)(s2*1e5) + (int)(s3*1e5/c) / 1e5 +
         graph.nNeighbors(i) / 1e9;
      q[qsize++] = make_pair(-weight, i);
   }
   sort(q, q + qsize);

   vector<int> ret;
   for (int i = 0; i < size && i < qsize; i++) {
      ret.push_back(q[i].second);
   }

   delete[] tp;
   delete[] heap;
   delete[] pos;
   for (unsigned i = 0; i < graph.nNodes; i++) {
      delete[] sh[i];
      delete[] influential[i];
   }
   delete[] sh;
   delete[] influential;

   return ret;
}

vector<int>
HIS(int k, double precision, vector<int> target_communities, vector<int> n2c,
    double * page_rank, const Graph & graph)
{
   int c         = target_communities.size();
   double alpha  = 0.3;
   double * beta = new double[two(c)];

   for (int set = 0; set < two(c); set++) {
      if (countbit(set)==0) {
         beta[set]=0;
      } else if (countbit(set)==1) {
         beta[set]=0;
      } else if (countbit(set)==2) {
         beta[set]=0.17;
      } else if (countbit(set)==3) {
         beta[set]=0.25;
      } else if (countbit(set)==4) {
         beta[set]=0.29;
      } else if (countbit(set)==5) {
         beta[set]=0.30;
      } else if (countbit(set)>=6) {
         beta[set]=0.35;
      }
   }

   vector<int> r2 = structure_hole_min_max_faster(c, k, precision, alpha, beta,
                                                  page_rank, target_communities,
                                                  n2c, graph);

   delete[] beta;
   
   return r2;
}

map<int, vector<vector<int> > >
allCliquesRec(const Graph & graph, int currentNode, vector<int> & nodes,
              int minSize, int maxSize)
{
   map<int, vector<vector<int> > > cliques;

   if (nodes.size() > maxSize)
      return cliques;

   if(currentNode == graph.nNodes) {
      if (nodes.size() >= minSize) {
         bool isClique = true;
         for(int i = 0; i < nodes.size(); i++){
            for(int j = i + 1; j < nodes.size(); j++){
               auto p = graph.neighbors(nodes[i]);
               bool areNeighbors = false;
               for (int k = 0; k < graph.nNeighbors(nodes[i]); k++) {
                  if (*(p.first + k) == nodes[j]) {
                     areNeighbors = true;
                  }
               }
               if(areNeighbors == false) {
                  isClique = false;
               }
            }
         }
         if(isClique == true){
            cliques[nodes.size()].push_back(nodes);
         }
      }
      return cliques;
   }

   //if x < n
   map<int , vector<vector<int> > > tmp;
   tmp = allCliquesRec(graph, currentNode + 1, nodes, minSize, maxSize);
   for (int i = minSize; i <= maxSize; i++) {
      cliques[i].insert(cliques[i].end(), tmp[i].begin(), tmp[i].end());
   }

   nodes.push_back(currentNode);

   tmp = allCliquesRec(graph, currentNode + 1, nodes, minSize, maxSize);
   for (int i = minSize; i <= maxSize; i++) {
      cliques[i].insert(cliques[i].end(), tmp[i].begin(), tmp[i].end());
   }

   nodes.pop_back();

   return cliques;
}

map<int, vector<vector<int> > >
allCliques(const Graph & graph, int minSize, int maxSize)
{
   vector<int> nodes;
   return allCliquesRec(graph, 0, nodes, minSize, maxSize);
}

vector<float>
weightedKCore(const Graph & graph, int nCores, vector<int> & outNodeToCore)
{
   outNodeToCore.resize(graph.nNodes, -1);

   float maxDegree = 0;
   vector<float> degrees(graph.nNodes, 0.);

   for (int idNode = 0; idNode < graph.nNodes; idNode++) {
      degrees[idNode] = graph.weightedDegree(idNode);

      if (degrees[idNode] > maxDegree) {
         maxDegree = degrees[idNode];
      }
   }

   vector<int> nodes;
   for (int idNode = 0; idNode < graph.nNodes; idNode++) {
      nodes.push_back(idNode);
   }
   sort(nodes.begin(), nodes.end(), [&](int i1, int i2){ return degrees[i1] < degrees[i2]; });

   int currentBinSize = 0;
   int idCore         = 0;
   int sizeBin        = (int)(graph.nNodes / nCores);

   for (int currentNodeId = 0; currentNodeId < graph.nNodes; currentNodeId++) {
      int node = nodes[currentNodeId];

      currentBinSize++;
      if (currentBinSize >= sizeBin && idCore < nCores - 1) {
         currentBinSize = 0;
         idCore++;
      } 
      outNodeToCore[node] = idCore;

      auto p = graph.neighbors(node);
      for (int idNgb = 0; idNgb < graph.nNeighbors(node); idNgb++) {
         unsigned other = *(p.first + idNgb);
         unsigned cost  = *(p.second + idNgb);

         if (degrees[other] > degrees[node]) {
            degrees[other] -= cost;

            int idxNode     = distance(nodes.begin(), find(nodes.begin() + currentNodeId, nodes.end(), other));
            int idxShouldBe = idxNode;

            while (idxShouldBe < graph.nNodes &&
                   degrees[nodes[idxNode]] > degrees[nodes[idxShouldBe]]) {
               idxShouldBe++;
            }
            
            if (idxNode != idxShouldBe) {
               int tmp = nodes[idxNode];
               nodes.erase(nodes.begin() + idxNode);
               nodes.insert(nodes.begin() + idxShouldBe -1, tmp);
            }
         }
      }
   }

   return degrees;
}
