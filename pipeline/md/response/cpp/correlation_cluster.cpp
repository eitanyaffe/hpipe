#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <math.h>

#include <queue>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "Params.h"

using namespace std;

enum Method { mSpearman, mPearson };

void init_params(int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("input table"), true);
  params.add_parser("ofn", new ParserFilename("output table"), true);
  params.add_parser("method", new ParserString("correlation function"), true);
  params.add_parser("ofn_bins", new ParserFilename("output pearson bin table"), true);
  params.add_parser("threshold", new ParserDouble("correlation linkage threshold", 0.95), false);
  params.add_parser("min_score", new ParserDouble("minimal coverage skew score", -4), false);

  if (argc == 1) {
    params.usage(argv[0]);
    exit(0);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I/O
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double read_double(string fn)
{
  cout << "reading file: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  string line;
  getline(in, line);
  double result = atof(line.c_str());
  in.close();
  return result;
}

void read_table(string fn,
		vector< vector<double> >& table, bool to_rank, double min_score,
		map<int, string>& map_back,
		int& N, int& M)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  M = fields.size() - 1;

  N = 0;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string contig = fields[0];
    // compute x_j - mean(x)
    vector<double> vals(M);
    double sum = 0;
    for (int j=0; j<M; ++j) {
      vals[j] = atof(fields[j+1].c_str());
      if (vals[j] < min_score)
	vals[j] = min_score;
      sum += vals[j];
    }
    double mean = sum / (double)M;
    for (int j=0; j<M; ++j)
      vals[j] -= mean;

    if (to_rank) {
      vector<double> sorted = vals;
      sort(sorted.begin(), sorted.end());
      map<double, int> rank_map;
      for (int i=0; i<M; i++) {
	double v = sorted[i];
	if (rank_map.find(v) == rank_map.end())
	  rank_map[v] = i;
      }

      vector<double> ranks(M);
      for (int i=0; i<M; i++) {
	massert(rank_map.find(vals[i]) != rank_map.end(), "value %f not found in map", vals[i]);
	ranks[i] = rank_map[vals[i]];
	// cout << "i=" << i << ", v=" << vals[i] << ", rank=" << ranks[i] << endl;
      }
      // exit(1);
      vals = ranks;
    }
    table.push_back(vals);
    map_back[N++] = contig;
  }
}

inline double get_pearson(const vector< vector<double> >& table, int i, int j, int M)
{
  double a_i=0, a_j=0, a_both=0;
  for (int m=0; m<M; m++) {
    a_i += table[i][m] * table[i][m];
    a_j += table[j][m] * table[j][m];
    a_both += table[i][m] * table[j][m];
  }
  double div = sqrt(a_i) * sqrt(a_j);
  if (div == 0)
    return (-2);
  else
    return (a_both / div);
}

inline double get_spearman(const vector< vector<double> >& table, int i, int j, int M)
{
  double r = 0;
  for (int m=0; m<M; m++) {
    double d = table[i][m] - table[j][m];
    r += d*d;
  }
  return (1 - (6 * r) / (M * (M*M-1)));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  Parameters params;
  init_params(argc, argv, params);

  string ifn = params.get_string("ifn");
  double threshold = params.get_double("threshold");
  double min_score = params.get_double("min_score");
  string ofn = params.get_string("ofn");
  string ofn_bins = params.get_string("ofn_bins");
  string method_str = params.get_string("method");
  massert(method_str == "pearson" || method_str == "spearman", "supported methods are 'pearson' or 'spearman'");
  Method method = method_str == "pearson" ? mPearson : mSpearman;

  vector< vector<double> > table; // contig[condition[i]]
  map<int, string> map_back;      // map from contig index back to name for output
  int N;                        // number of contigs
  int M;                        // number of conditions per contig

  read_table(ifn, table, method == mSpearman, min_score, map_back, N, M);

  cout << "number of contigs: " << N << endl;
  cout << "number of conditions per contig: " << M << endl;

  vector< set<int> > clusters;
  map < int, int> contig2cluster;

  // keep track of removed clusters
  map< int, bool > removed;

  // bin for hist
  int bin_count = 200;
  vector<double> bins(bin_count);
  cout << "number of bins: " << bins.size() << endl;
  for (unsigned int i=0; i<bins.size(); i++)
    bins[i] = 0;

  for (int i=0; i<N; i++) {
    if (i > 0 && i % 10000 == 0)
      cout << "contig progress: " << i << endl;

    // compute all neighbours of contig_i among earlier contigs
    set<int> neighbours;
    for (int j=0; j<i; j++) {
      double rho = ((method == mPearson) ? get_pearson(table, i, j, M) : get_spearman(table, i, j, M));
      if (rho == -2)
	continue;
      if ((map_back[i] == "c653965" && map_back[j] == "c9527") || (map_back[j] == "c653965" && map_back[i] == "c9527"))
	cout << "rho=" << rho << endl;

      unsigned int bin = (int)((rho/2+0.5) * bin_count);
      if (bin < 0) bin = 0;
      if (bin > bins.size()-1) bin = bins.size()-1;
      bins[bin]++;

      // cout << "rho: " << rho << endl;
      if (rho > threshold) {
	massert(contig2cluster.find(j) != contig2cluster.end(), "contig not found: %d", j);
	unsigned int cluster_index = contig2cluster[j];
	neighbours.insert(cluster_index);
      }
    }
    // cout << "n_clusters=" << clusters.size() << ", n_neighbours=" << neighbours.size() << endl;

    if (neighbours.size() == 0) {
      // new singlton
      set<int> cluster;
      cluster.insert(i);
      int cluster_index = clusters.size();
      contig2cluster[i] = cluster_index;
      clusters.push_back(cluster);
    } else if (neighbours.size() == 1) {
      // add to single cluster found
      unsigned int cluster_index = (*neighbours.begin());
      massert(cluster_index >= 0 && cluster_index < clusters.size(), "cluster out of range: %d", cluster_index);
      set<int>& cluster = clusters[cluster_index];
      cluster.insert(i);
      contig2cluster[i] = cluster_index;
    } else {
      // unite all previous clusters into new cluster
      int new_cluster_index = clusters.size();
      set<int> new_cluster;

      for (set<int>::iterator it = neighbours.begin(); it != neighbours.end(); ++it) {
	unsigned int ocluster_index = *it;
	massert(ocluster_index >= 0 && ocluster_index < clusters.size(), "cluster out of range: %d", ocluster_index);
	set<int>& ocluster = clusters[ocluster_index];
	for (set<int>::iterator jt = ocluster.begin(); jt != ocluster.end(); ++jt) {
	  int oocontig = (*jt);
	  new_cluster.insert(oocontig);
	  contig2cluster[oocontig] = new_cluster_index;
	}
	// invalidate old cluster
	removed[ocluster_index] = true;
      }

      // finally add the new contig
      new_cluster.insert(i);
      contig2cluster[i] = new_cluster_index;

      clusters.push_back(new_cluster);
    }
  }

  cout << "total clusters generated (including invalidated clusters): " << clusters.size() << endl;
  map<int, int> sizes;
  for (int i=0; i<N; i++) {
    massert(contig2cluster.find(i) != contig2cluster.end(), "contig not found: %d", i);
    int cluster_idx = contig2cluster[i];
    sizes[cluster_idx]++;
  }

  // output result
  cout << "writing output cluster file: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  map<int,int> final_index_map;
  out << "contig\tcluster\tinternal_index" << endl;
  int new_index = 1;
  for (int i=0; i<N; i++) {
    massert(contig2cluster.find(i) != contig2cluster.end(), "contig not found: %d", i);
    int cluster_idx = contig2cluster[i];
    int size = sizes[cluster_idx];
    if (size == 1) {
      out << map_back[i] << "\t" << -1 << endl;
    } else {
      if (final_index_map.find(cluster_idx) == final_index_map.end())
	final_index_map[cluster_idx] = new_index++;
      int index = final_index_map[cluster_idx];
      out << map_back[i] << "\t" << index << "\t" << cluster_idx << endl;
    }
  }
  out.close();

  int scount = 0;
  int lcount = 0;
  for(map<int, int>::iterator it = sizes.begin(); it != sizes.end(); ++it) {
    int size = (*it).second;
    if (size == 1)
      scount++;
    else
      lcount++;
  }
  cout << "singlton count: " << scount << endl;
  cout << ">1 count: " << lcount << endl;

  // output bins
  cout << "writing bins: " << ofn_bins << endl;
  ofstream out2(ofn_bins.c_str());
  massert(out2.is_open(), "could not open file %s", ofn_bins.c_str());
  out2 << "bin\tcount" << endl;
  for (unsigned int i=0; i<bins.size(); i++)
    out2 << (double)i/(bin_count/2) - 1 << "\t" << bins[i] << endl;
  out.close();

}
