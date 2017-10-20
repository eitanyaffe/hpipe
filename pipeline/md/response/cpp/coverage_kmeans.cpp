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
  params.add_parser("ofn_prefix", new ParserFilename("output prefix"), true);
  params.add_parser("n_clusters", new ParserInteger("number of clusters"), true);
  params.add_parser("random_seed", new ParserInteger("random generator seed", 0), false);
  params.add_parser("cluster_seed_method", new ParserString("cluster seed method [diameter|random]", "random"), false);

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

void read_table(string fn,
		vector< vector<double> >& table,
		map<int, string>& map_back,
		unsigned int& N, unsigned int& M)
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

    // contig
    string contig = fields[0];

    // coverage
    vector<double> vals(M);
    double total = 0;
    for (unsigned int j=0; j<M; ++j) {
      vals[j] = atof(fields[j+1].c_str());
      total += vals[j];
    }
    table.push_back(vals);

    // map back to contig string identifier
    map_back[N++] = contig;
  }
}

inline double get_pearson(const vector<double>& v1, const vector<double>& v2)
{
  massert(v1.size() == v2.size(), "vectors must be same size");
  int N = v1.size();

  double a_1=0, a_2=0, a_both=0;
  for (int i=0; i<N; i++) {
    a_1 += v1[i] * v1[i];
    a_2 += v2[i] * v2[i];
    a_both += v1[i] * v2[i];
  }
  double div = sqrt(a_1) * sqrt(a_2);
  if (div == 0)
    return (-2);
  else
    return (a_both / div);
}

inline double distance(const vector<double>& v1, const vector<double>& v2)
{
  massert(v1.size() == v2.size(), "vectors must be same size");
  int N = v1.size();

  double a=0;
  for (int i=0; i<N; i++)
    a += (v1[i] - v2[i]) * (v1[i] - v2[i]);

  return sqrt(a);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  Parameters params;
  init_params(argc, argv, params);

  string ifn = params.get_string("ifn");
  unsigned int n_clusters = params.get_int("n_clusters");
  string ofn_prefix = params.get_string("ofn_prefix");
  unsigned int rnd_seed = params.get_int("random_seed");
  string cluster_seed_method = params.get_string("cluster_seed_method");

  // set random seed
  srand(rnd_seed);

  vector< vector<double> > table; // contig[condition[i]]
  map<int, string> map_back;      // map from contig index back to name for output
  unsigned int n_contigs;         // number of contigs
  unsigned int n_conditions;      // number of conditions per contig
  read_table(ifn, table, map_back, n_contigs, n_conditions);

  cout << "number of clusters: " << n_clusters << endl;
  cout << "number of contigs: " << n_contigs << endl;
  cout << "number of conditions per contig: " << n_conditions << endl;

  vector< vector<int> > cluster_contigs(n_clusters);
  vector<bool> is_dead(n_clusters);
  vector< vector<double> > centroids(n_clusters);

  for (unsigned int j=0; j<n_clusters; j++)
     is_dead[j] = false;

  ////////////////////////////////////////////////////////////////////////////////////////////
  // set seed centroids
  ////////////////////////////////////////////////////////////////////////////////////////////

  for (unsigned int j=0; j<n_clusters; j++) {
    int contig_index;

    // get farthest contig index
    if (j == 0 || cluster_seed_method == "random") {
      contig_index = rand() % n_contigs;
    } else {
      double max_dist = -1;
      for (unsigned int i=0; i<n_contigs; i++) {

	// find minimal distance to all previous defined seed centroids
	double cent_dist = -1;
	for (unsigned int k=0; k<j; k++) {
	  double dist = distance(centroids[k], table[i]);
	  if (dist < cent_dist || cent_dist == -1) {
	    cent_dist = dist;
	  }
	}

	if (cent_dist > max_dist || i == 0) {
	  max_dist = cent_dist;
	  contig_index = i;
	}
      }

      // !!!
      cout << "max_dist: " << max_dist << endl;
    }


    // set centroid
    centroids[j] = table[contig_index];


    // print seed centroid
    cout << "seed " << j << ", contig " << map_back[contig_index] << " : ";
    for (unsigned int m=0; m<n_conditions; m++)
      cout << ", " << centroids[j][m];
    cout << endl;
  }

  vector<int> clusters(n_contigs);
  for (unsigned int i=0; i<n_contigs; i++)
    clusters[i] = -1;


  ////////////////////////////////////////////////////////////////////////////////////////////
  // k-means core
  ////////////////////////////////////////////////////////////////////////////////////////////

  int iteration = 1;
  int max_iteration = 10000;

  while (1) {

    // how many contigs changed cluster?
    int n_changed = 0;

    // erase previous cluster_contigs
    for (unsigned int j=0; j<n_clusters; j++)
      cluster_contigs[j].clear();

    // 1: assignment step
    for (unsigned int i=0; i<n_contigs; i++) {
      int closest_cluster = -1;
      int min_dist = -1;

      // !!!
      // cout << "i: " << i << endl;

      for (unsigned int j=0; j<n_clusters; j++) {
	if (is_dead[j])
	  continue;

	double dist = distance(centroids[j], table[i]);

	// !!!
	// cout << "j: " << j << ", dist=" << dist << endl;

	if (closest_cluster == -1 || dist < min_dist) {
	  closest_cluster = j;
	  min_dist = dist;
	}
      }
      massert(closest_cluster != -1, "no cluster found");
      n_changed += (clusters[i] != closest_cluster) ? 1 : 0;

      clusters[i] = closest_cluster;
      cluster_contigs[closest_cluster].push_back(i);

      // !!!
      // cout << "closest: " << closest_cluster << endl;
      // if (i > 10)
      // 	exit(1);
    }

    cout << "iteration: " << iteration++ << ", n_changed: " << n_changed << endl;

    if (iteration > max_iteration) {
      cout << "Warning: stability not reached" << endl;
      break;
    }

    // stop when stable
    if (n_changed == 0)
      break;

    // 2: update step
    for (unsigned int j=0; j<n_clusters; j++) {
      unsigned int n_cluster_size = cluster_contigs[j].size();
      // cout << "cluster " << j << " size=" << n_cluster_size << endl;
      is_dead[j] = (n_cluster_size == 0);
      if (is_dead[j])
	continue;

      vector<double> new_centroid(n_conditions);

      // init empty centroid
      for (unsigned int m=0; m<n_conditions; m++)
	new_centroid[m] = 0;

      // sum up all cluster members
      for (unsigned int i=0; i<n_cluster_size; i++) {
	for (unsigned int m=0; m<n_conditions; m++) {
	  int contig_index = cluster_contigs[j][i];
	  new_centroid[m] += table[contig_index][m];
	}
      }

      // divide to get mean
      for (unsigned int m=0; m<n_conditions; m++)
	new_centroid[m] /= n_cluster_size;

      centroids[j] = new_centroid;
    }

  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  // output
  ////////////////////////////////////////////////////////////////////////////////////////////

  string ofn_clusters = ofn_prefix + ".clusters";
  string ofn_centroid_mean = ofn_prefix + ".centroid_mean";
  string ofn_centroid_sd = ofn_prefix + ".centroid_sd";

  // output cluster assignments
  cout << "writing clusters to file: " << ofn_clusters << endl;
  ofstream out(ofn_clusters.c_str());
  massert(out.is_open(), "could not open file %s", ofn_clusters.c_str());
  out << "contig\tcluster" << endl;
  for (unsigned int i=0; i<n_contigs; i++)
    out << map_back[i] << "\t" << clusters[i] << endl;
  out.close();

  // output centroid mean
  cout << "writing centroid means to file: " << ofn_centroid_mean << endl;
  ofstream out2(ofn_centroid_mean.c_str());
  massert(out2.is_open(), "could not open file %s", ofn_centroid_mean.c_str());
  out2 << "centroid\tsize";
  for (unsigned int m=0; m<n_conditions; m++)
    out2 << "\t"<< "condition_" << m;
  out2 << endl;
  for (unsigned int j=0; j<n_clusters; j++) {
    out2 << j << "\t"  << cluster_contigs[j].size();
    for (unsigned int m=0; m<n_conditions; m++)
      out2 << "\t" << centroids[j][m];
    out2 << endl;
  }
  out2.close();

  // output centroid sd
  cout << "writing centroid sds to file: " << ofn_centroid_sd << endl;
  ofstream out3(ofn_centroid_sd.c_str());
  massert(out3.is_open(), "could not open file %s", ofn_centroid_sd.c_str());
  out3 << "centroid\tsize";
  for (unsigned int m=0; m<n_conditions; m++)
    out3 << "\t"<< "condition_" << m;
  out3 << endl;
  for (unsigned int j=0; j<n_clusters; j++) {
    unsigned int n_cluster_size = cluster_contigs[j].size();
    out3 << j << "\t"  << n_cluster_size;
    for (unsigned int m=0; m<n_conditions; m++) {
      double total = 0;
      for (unsigned int i=0; i<n_cluster_size; i++) {
	int contig_index = cluster_contigs[j][i];
	double v = centroids[j][m] - table[contig_index][m];
	total += v * v;
      }
      double sd = n_cluster_size > 0 ? sqrt(total / n_cluster_size) : 0;
      out3 << "\t" << sd;
    }
    out3 << endl;
  }
  out3.close();
}
