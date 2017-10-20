#ifndef __MODEL__
#define __MODEL__

#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <math.h>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

// Access matrix in vector form
#define MATRIX_INDEX(x, y, dim) ((x-1) + (y-1)*dim)

using namespace std;

// The model has a set of features. Each feature is factor which is a function of
// the two interacting fends
struct ModelFeatures
{
  int size;
  vector<string> fields;
  vector<string> fns;
  vector<int> sizes;
};

struct Fend
{
  int index;
  int contig_index;
  int coord;
  vector<int> levels;
};

struct comp_fends : public binary_function<Fend, Fend, bool> {
  bool operator()(Fend x, Fend y) { return x.index < y.index; }
};

struct Contig
{
  string id; // contig id
  int anchor; // index of anchor, -1 of none
  Contig(string _id, int _anchor) : id(_id), anchor(_anchor) {};
};

class Model
{
 private:

  // keep fends by contig index and anchor
  map<int, vector<Fend> > m_fends_contig;
  map<int, vector<Fend> > m_fends_anchor;


  // features
  ModelFeatures m_features;

  // prior
  double m_prior;

  // correction matrices
  vector < vector<double> > m_model_matrices;

  // map from fend to contig id
  map<int, string> m_fend_contig_map;

  // map from contig to anchor
  //  map<string, int> m_contig_anchor_map;

  vector<Contig> m_contigs;
  map<string, bool> m_contigs_flag;

  // map from contig id to contig index and back
  map<string, int> m_contig_index_map;
  map<int, string> m_rev_contig_index_map;
  int contig2int(string contig) {
    if (m_contig_index_map.find(contig) == m_contig_index_map.end())
      {
	int index = m_contig_index_map.size();
	m_contig_index_map[contig] = index;
	m_rev_contig_index_map[index] = contig;
      }
    return m_contig_index_map[contig];
  };

  void read_feature_table(string fn, vector<double>& result, int mat_dim);
  void read_feature_tables(const vector<string>& feature_fns,
			   const vector<int>& sizes);

  void intersect_fends(vector<Fend>& fends_x, vector<Fend>& fends_y, vector<Fend>& fends_common);
  double compute_simple(vector<Fend>& fends_x, vector<Fend>& fends_y);

 public:

  Model(ModelFeatures features, double prior);

  void read_fends_table(const string& fn,
			const string& fend_field,
			const string& contig_field,
			const string& cooord_field,
			const string& anchor_field,
			int& n_anchors);

  // map from fend index to contig id
  map<int, string>& fend_contig_map() { return m_fend_contig_map; };

  // map from contig id to anchor
  //  map<string, int>& contig_anchor_map() { return m_contig_anchor_map; };

  vector<Contig>& contigs() { return m_contigs; };

  map<string, int>& contig_index_map() { return m_contig_index_map; };
  map<int, string>& rev_contig_index_map() { return m_rev_contig_index_map; };

  // compute expected contacts between pair of contigs
  double compute(vector<Fend>& fends_x, vector<Fend>& fends_y);
  double compute(string contig, int anchor);
  double compute(string contig1, string contig2);
};

#endif
