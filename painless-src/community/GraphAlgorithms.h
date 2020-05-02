#pragma once

#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

#include "Graph.h"

using namespace std;

#define two(X) (1<<(X))
#define contain(S,X) (((S)&two(X))!=0)

static const double pi  = acos(-1.0);

template<class T>
static inline void checkmax(T &a,T b)
{
   if(b > a) a = b;
}

vector<map<int, double> > renumber(const vector< map<int, double> > & graph,
                                    vector<int> & r2o, vector<int> & o2r);

void pageRank(const Graph & graph, double * page_rank);

void weightedPageRank(const Graph & graph, double * page_rank);

vector<int> HIS(int k, double precision, vector<int> target_communities,
                vector<int> n2c, double * page_rank, const Graph & graph);

void louvain(const Graph & graph, double precision, const vector<int> & r2o,
             vector<vector<int> > & area2nodes, vector<int> & node2area);

map<int, vector<vector<int> > > allCliques(const Graph & graph, int minSize,
                                           int maxSize);

vector<float> weightedKCore(const Graph & graph, int nCores,
                            vector<int> & outNodeToCore);
