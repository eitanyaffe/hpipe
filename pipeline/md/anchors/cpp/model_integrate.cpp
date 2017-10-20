#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <math.h>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

// Access matrix in vector form
#define MATRIX_INDEX(x, y, dim) ((x-1) + (y-1)*dim)

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/*

  Input
  -----
  1. A table of fragment ends (fends). Each fend has a list of features
  (such as frag_len_bin, gc_bin, etc.)

  2. A list of 0 or more inter/intra anchor model features. Each model feature is defined by:
  - field name: same name as the one in the fend table
  - file: contains correction matrix for feature
  - size: matrix dimension

  3. Two list (for dimension x and y) of 1 or more output features. Each feature is defined by:
  - field name: same name as the one in the fend table
  - from: feature start index
  - to: feature to index

  4. An inter-anchor prior

  5. An intra-anchor prior table

  Output
  ------
  The program outputs a single file that contains model-expected counts
  for the requested output bins.

*/
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures
//////////////////////////////////////////////////////////////////////////////////////////////////

// The model has a set of features. Each feature is factor which is a function of
// the two interacting fends
struct ModelFeatures
{
  int size;
  vector<string> fields;
  vector<int> sizes;

  vector<string> fns;
};

// We integrate over two sets of bins, one for each dimension (x,y)
struct OutputFeatures
{
  int size;
  vector<string> fields;
  vector<int> from;
  vector<int> to;
};

// every fend is defined by a set of levels, one for each model feature
struct Fend
{
  int index;
  int contig_code;
  int anchor;
  bool on_anchor;

  vector<int> levels;
};

// every bin holds a collection of all fends that have the same levels for all model features
struct Bin
{
  vector<Fend> fends;
};

// limit fend pairs to score
enum Scope { s_anchored, s_inter_anchor, s_intra_anchor };
Scope string2scope(string str)
{
  Scope scope;
  if (str == "inter_anchor")
    scope = s_inter_anchor;
  else if (str == "intra_anchor")
    scope = s_intra_anchor;
  else if (str == "anchored")
    scope = s_anchored;
  else {
    cerr << "Unknown scope " << str << endl;
    exit(1);
  }
  return (scope);
}
string scope2string(Scope scope)
{
  string result;
  switch (scope) {
  case s_inter_anchor:
    result = "inter_anchor";
    break;
  case s_intra_anchor:
    result = "intra_anchor";
    break;
  case s_anchored:
    result = "anchored";
    break;
  default:
    cerr << "Unknown scope index " << scope << endl;
    exit(1);
  }
  return (result);
}

// model
enum Model { m_count, m_std };
Model string2model(string str)
{
  Model model;
  if (str == "count")
    model = m_count;
  else if (str == "std")
    model = m_std;
  else {
    cerr << "Unknown model " << str << endl;
    exit(1);
  }
  return (model);
}
string model2string(Model model)
{
  string result;
  switch (model) {
  case m_count:
    result = "count";
    break;
  case m_std:
    result = "std";
    break;
  default:
    cerr << "Unknown model index " << model << endl;
    exit(1);
  }
  return (result);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Parsing and handling user arguments
//////////////////////////////////////////////////////////////////////////////////////////////////

class UserParams {
public:
  string fends_fn;
  string ofn;
  // limit output to fend range
  int from_fend, to_fend;

  // scope of fend pairs
  Scope scope;

  // model
  Model model;

  // model parameters
  double prior;
  ModelFeatures m_features;

  // output features
  OutputFeatures o_features_x;
  OutputFeatures o_features_y;

  // ctor with default values
  UserParams() {
    fends_fn = "results/wtb_s0_trans.binned";
    ofn = "out";
    from_fend = 0;
    to_fend = 0;
    scope = s_anchored;
    model = m_count;
    prior = 1;
  }
};

ostream& operator <<(ostream &os, const UserParams &params)
{
  os << "binned fends table: " << params.fends_fn << endl;
  os << "output file: " << params.ofn << endl;

  os << "fend from: " << params.from_fend << endl;
  os << "fend to: " << params.to_fend << endl;

  os << "scope: " << scope2string(params.scope) << endl;
  os << "model: " << model2string(params.model) << endl;
  os << "prior: " << params.prior << endl;

  // model features
  if (params.m_features.size == 0)
    os << "No model functions" << endl;
  else for (int i=0; i<params.m_features.size; i++)
	 os << "M" << i << ": " << params.m_features.fields[i]
	    << ", dim=" << params.m_features.sizes[i] << endl;

  // output features
  for (int i=0; i<params.o_features_x.size; i++)
    cerr << "Output bin X" << i << ": " << params.o_features_x.fields[i] << ", " <<
      params.o_features_x.from[i] << "-" << params.o_features_x.to[i] << endl;
  for (int i=0; i<params.o_features_y.size; i++)
    cerr << "Output bin Y" << i << ": " << params.o_features_y.fields[i] << ", " <<
      params.o_features_y.from[i] << "-" << params.o_features_y.to[i] << endl;

  os << "-----------------------------------" << endl;

  return os;
}

void usage(const char* name, UserParams& params)
{
  fprintf(stderr, "usage: %s [options]\n", name);
  cerr << " -fends <fn>: input binned fend table" << endl;
  cerr << " -from_fend <fend index>: fend mask start" << endl;
  cerr << " -to_fend <fend index>: fend mask end" << endl;
  cerr << " -scope <str>: limit fend pairs to anchored|intra_anchor|inter_anchor contacts" << endl;
  cerr << " -model <str>: type of model count|std" << endl;
  cerr << " -model_num <n>: number of model features" << endl;
  cerr << " -f_model_field <string>: model feature field" << endl;
  cerr << " -f_model_size <n>: model feature matrix dimension" << endl;
  cerr << " -f_model_fn <filename>: model feature filename" << endl;
  cerr << " -output_{x|y}_num <n>: number of output x|y features" << endl;
  cerr << " -f_{x|y}_field <string>: output x|y feature field" << endl;
  cerr << " -f_{x|y}_from <int>: output x|y from bin index" << endl;
  cerr << " -f_{x|y}_to <int>: output x|y to bin index" << endl;
  cerr << " -output <fn>: output file" << endl;

  exit(1);
}

void parse_user_arguments(int argc, char **argv, UserParams& params)
{
  if (argc == 1)
    usage(argv[0], params);

  int i = 1;
  while (i < argc)
    {
      string option = argv[i];
      char* arg = argv[i+1];

      if (option == "-fends")
	params.fends_fn = arg;
      else if (option == "-output")
	params.ofn = arg;
      else if (option == "-from_fend")
	params.from_fend = atoi(arg);
      else if (option == "-to_fend")
	params.to_fend = atoi(arg);
      else if (option == "-scope")
	params.scope = string2scope(string(arg));
      else if (option == "-model")
	params.model = string2model(string(arg));
      else if (option == "-prior")
	params.prior = atof(arg);
      else if (option == "-f_model_field")
	params.m_features.fields.push_back(arg);
      else if (option == "-f_model_size")
	params.m_features.sizes.push_back(atoi(arg));
      else if (option == "-f_model_fn")
	params.m_features.fns.push_back(arg);
      else if (option == "-model_num")
	params.m_features.size = atoi(arg);
      else if (option == "-f_x_field")
	params.o_features_x.fields.push_back(arg);
      else if (option == "-f_x_from")
	params.o_features_x.from.push_back(atoi(arg));
      else if (option == "-f_x_to")
	params.o_features_x.to.push_back(atoi(arg));
      else if (option == "-output_x_num")
	params.o_features_x.size = atoi(arg);
      else if (option == "-f_y_field")
	params.o_features_y.fields.push_back(arg);
      else if (option == "-f_y_from")
	params.o_features_y.from.push_back(atoi(arg));
      else if (option == "-f_y_to")
	params.o_features_y.to.push_back(atoi(arg));
      else if (option == "-output_y_num")
	params.o_features_y.size = atoi(arg);
      else {
	cerr << "Error: unknown option: " << option << endl;
	exit(1);
      }

      i += 2;
    }

  // some sanity checks
  if (params.m_features.size != (int)params.m_features.fields.size()) {
    cerr << "Error: Expecting " << params.m_features.size << " model feature fields" << endl;
    exit(1);
  }
  if (params.m_features.size != (int)params.m_features.sizes.size()) {
    cerr << "Error: Expecting " << params.m_features.size << " model feature sizes" << endl;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Low level utility functions
//////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
string join_strings(T it, T end, string delim, string suffix)
{
  string result = *it + suffix;

  while (++it != end)
    result += delim + (*it) + suffix;
  return result;
}

template <typename T>
string NumberToString ( T Number )
{
  stringstream ss;
  ss << Number;
  return ss.str();
}

void massert(bool cond, char *fmt, ...)
{
  if (cond)
    return;

  fprintf(stderr, "Error: ");

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  fprintf(stderr, "\n");
  exit(-1);
}

// use a global structure for contig codes
static map<string, int> contig_codes;
static map<int, string> rev_contig_codes;
static int next_code = 0;
int contig2code(string chr)
{
  if (contig_codes.find(chr) == contig_codes.end())
    {
      contig_codes[chr] = next_code;
      rev_contig_codes[next_code] = chr;
      next_code++;
    }
  return contig_codes[chr];
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// I/O utility functions
//////////////////////////////////////////////////////////////////////////////////////////////////

void split_line(istream &in, vector<string> &fields, char delim)
{
  fields.resize(0);
  string field;
  while(in) {
    char c = in.get();
    if(c == '\r') {
      continue;
    }
    if(c == '\n') {
      fields.push_back(field);
      break;
    }
    if(c == delim) {
      fields.push_back(field);
      field.resize(0);
    } else {
      field.push_back(c);
    }
  }
}

string levels_to_string(const vector<int>& levels)
{
  assert(levels.size() > 0);
  string result = NumberToString(levels[0]);
  for (unsigned int i = 1; i<levels.size(); i++)
    result += string("\t") + NumberToString(levels[i]);
  return result;
}

void parse_levels(vector<int>& levels, const vector<int>& fields_ind, const vector<string>& fields)
{
  levels.resize(fields_ind.size());
  for (unsigned int i = 0; i<fields_ind.size(); i++)
    levels[i] = atoi(fields[fields_ind[i]].c_str());
}

void add_fend_to_obins(const Fend& fend,
                       map< string, Bin > &obins,
                       const vector<int>& obin_fields_ind,
                       const vector<string>& fields)
{
  vector<int> obin_levels;
  parse_levels(obin_levels, obin_fields_ind, fields);
  vector<Fend>& fends = obins[levels_to_string(obin_levels)].fends;
  fends.push_back(fend);
}

void init_field_indices(vector<int>& field_ind, const vector<string> &fields, const vector<string>& titles)
{
  field_ind.resize(fields.size(), -1);
  for (unsigned int i = 0; i<titles.size(); i++)
    for (unsigned int j = 0; j<fields.size(); j++)
      if (titles[i] == fields[j])
	field_ind[j] = i;
  for (unsigned int i = 0; i<field_ind.size(); i++)
    massert(field_ind[i] != -1, "field %s not found", fields[i].c_str());
}

int get_field_index(string field, const vector<string>& titles, const string& fn)
{
  int result = -1;
  for (unsigned int i = 0; i<titles.size(); i++)
    if (titles[i] == field)
      result = i;
  if (result == -1) {
    cerr << "could not find field " << field << " in table " << fn << endl;
  }
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read fend file
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_fends_table(const string& fn, const ModelFeatures& m_features,
                      const OutputFeatures& o_features_x, const OutputFeatures& o_features_y,
                      const int from_fend, const int to_fend,
                      map<string, Bin>& obins_x, map<string, Bin>& obins_y)
{
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  vector<int> m_features_ind(m_features.size, -1);
  vector<int> o_features_x_ind(o_features_x.size, -1);
  vector<int> o_features_y_ind(o_features_y.size, -1);

  // parse header
  split_line(in, fields, delim);
  init_field_indices(m_features_ind, m_features.fields, fields);
  init_field_indices(o_features_x_ind, o_features_x.fields, fields);
  init_field_indices(o_features_y_ind, o_features_y.fields, fields);

  int fend_ind = get_field_index("fend", fields, fn);
  int contig_ind = get_field_index("contig", fields, fn);
  int anchor_ind = get_field_index("anchor", fields, fn);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      return;

    Fend fend;
    fend.index = atoi(fields[fend_ind].c_str());
    fend.contig_code = contig2code(fields[contig_ind]);
    fend.anchor = atoi(fields[anchor_ind].c_str());
    fend.on_anchor = fend.anchor != 0;

    // classify fend according to its model levels
    parse_levels(fend.levels, m_features_ind, fields);

    // add fend to output bins according to the output levels
    add_fend_to_obins(fend, obins_x, o_features_x_ind, fields);

    if ( (from_fend == 0 && to_fend == 0) || ( fend.index >= from_fend && fend.index <= to_fend) )
      add_fend_to_obins(fend, obins_y, o_features_y_ind, fields);
  }

  in.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read prior table file
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_prior_table(string fn, map<int,double>& result)
{
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int anchor_ind = get_field_index("anchor", fields, fn);
  int prior_ind = get_field_index("prior", fields, fn);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      return;

    int anchor = atoi(fields[anchor_ind].c_str());
    double prior = atoi(fields[prior_ind].c_str());
    result[anchor] = prior;
  }
  in.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read feature file
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_feature_table(string fn, vector<double>& result, int mat_dim)
{
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';
  const unsigned int width = 3;
  bool first = true;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      return;
    if (first)
      {
	first = false;
	continue;
      }

    assert(fields.size() == width);

    int i = atoi(fields[0].c_str());
    int j = atoi(fields[1].c_str());
    double prob = atof(fields[2].c_str());

    int v1 = MATRIX_INDEX(i, j, mat_dim);
    int v2 = MATRIX_INDEX(j, i, mat_dim);

    result[v1] = prob;
    result[v2] = prob;
  }

  in.close();
}

void read_feature_tables(const vector<string>& feature_fns, const vector<int>& sizes,
			 vector < vector<double> >& features)
{
  massert(feature_fns.size() == sizes.size(), "model matrix filenames not supplied");
  for (unsigned int i=0; i < feature_fns.size(); i++)
    {
      vector<double>& feature = features[i];
      feature.resize(sizes[i]*sizes[i]);
      read_feature_table(feature_fns[i], feature, sizes[i]);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Debug
//////////////////////////////////////////////////////////////////////////////////////////////////

void dump_bins(map<string, Bin>& bins, string title)
{
  cerr << title << " bins" << endl;
  map<string, Bin>::iterator it = bins.begin();
  while(it != bins.end())
    {
      cerr << "Bin: " << (*it).first << endl;

      vector<Fend>& fends = (*it).second.fends;
      vector<Fend>::iterator jt = fends.begin();
      while(jt != fends.end())
        {
	  Fend& fend = *jt;
	  cerr << " " << fend.index << " " << fend.contig_code << " " << rev_contig_codes[fend.contig_code] << ": ";
	  for (unsigned int i = 0; i < fend.levels.size(); i++)
	    cerr << fend.levels[i] << " ";
	  cerr << endl;
	  jt++;
        }
      it++;
      cerr << "=============================================================" << endl;
    }
}

void dump_feature(const vector<double>& f, int mat_dim)
{
  cerr << "==============================================" << endl;
  cerr << "i\tj\tprob" << endl;
  for (unsigned int k = 0; k < f.size(); k++)
    {

      int i = k / mat_dim + 1;
      int j = k - (i-1) * mat_dim + 1;

      double prob = f[k];
      cerr << i << "\t" << j << "\t" << prob << endl;
    }
}

void dump_features(const vector < vector<double> >& features, const vector<int>& sizes)
{
  for (unsigned int i = 0; i < features.size(); i++)
    dump_feature(features[i], sizes[i]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Parse command line
//////////////////////////////////////////////////////////////////////////////////////////////////

char** parse_output_features(char **argv, int& no, OutputFeatures& obins)
{
  no = atoi(argv[0]);
  argv++;

  obins.size = no;
  obins.fields.resize(no);
  obins.from.resize(no);
  obins.to.resize(no);

  const int items_per_field = 3;
  for (int i = 0; i < no; i++)
    {
      obins.fields[i] = argv[i*items_per_field];
      obins.from[i] = atoi(argv[i*items_per_field+1]);
      obins.to[i] = atoi(argv[i*items_per_field+2]);
    }
  return (argv + items_per_field*no);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Fend list intersection
//////////////////////////////////////////////////////////////////////////////////////////////////

struct comp_fends : public binary_function<Fend, Fend, bool> {
  bool operator()(Fend x, Fend y) { return x.index < y.index; }
};

void intersect_fends(vector<Fend>& fends_x, vector<Fend>& fends_y,
                     vector<Fend>& fends_x_unique, vector<Fend>& fends_y_unique, vector<Fend>& fends_common)
{
  sort(fends_x.begin(), fends_x.end(), comp_fends());
  sort(fends_y.begin(), fends_y.end(), comp_fends());

  // X && Y
  set_intersection(fends_x.begin(), fends_x.end(),
		   fends_y.begin(), fends_y.end(), back_inserter(fends_common), comp_fends());

  // X \ Y
  set_difference(fends_x.begin(), fends_x.end(),
		 fends_y.begin(), fends_y.end(), back_inserter(fends_x_unique), comp_fends());

  // Y \ X
  set_difference(fends_y.begin(), fends_y.end(),
		 fends_x.begin(), fends_x.end(), back_inserter(fends_y_unique), comp_fends());
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

void init_keys(const OutputFeatures& features, int index, string key, vector<string>& keys)
{
  if (index == features.size) {
    keys.push_back(key);
    return;
  }
  for (int i = features.from[index]; i <= features.to[index]; i++)
    {
      string next_key = (index == 0) ? NumberToString(i) : key + string("\t") + NumberToString(i);
      init_keys(features, index+1, next_key, keys);
    }
}

double fends_expected(const vector<Fend>& fends_x, const vector<Fend>& fends_y,
                      Scope scope, Model model,
                      const vector < vector<double> >& model_matrices,
		      const vector<int>& msizes,
                      const double prior)
{
  double total = 0;
  vector<Fend>::const_iterator jt_x, jt_y;

  for(jt_x = fends_x.begin(); jt_x != fends_x.end(); jt_x++)
    for(jt_y = fends_y.begin(); jt_y != fends_y.end(); jt_y++)
      {
        const Fend& fend_x = *jt_x;
        const Fend& fend_y = *jt_y;

        // skip intra-contig contacts
        if (fend_x.index == fend_y.index || fend_x.contig_code == fend_y.contig_code)
	  continue;

        switch (scope) {
        case s_inter_anchor:
	  if (!fend_x.on_anchor || !fend_y.on_anchor || (fend_x.anchor == fend_y.anchor))
	    continue;
	  break;
        case s_intra_anchor:
	  if (!fend_x.on_anchor || !fend_y.on_anchor || (fend_x.anchor != fend_y.anchor))
	    continue;
	  break;
        case s_anchored:
	  if (!fend_x.on_anchor && !fend_y.on_anchor)
	    continue;
	  break;
        default:
	  cerr << "Error: unknown scope" << endl;
	  exit(1);
	}

	double value = 0;
        if (model == m_count) {
	  value = 1;
	} else if (model == m_std) {
	  value = prior;
	  for(unsigned int mi = 0; mi < model_matrices.size(); mi++)
	    value *= model_matrices[mi][MATRIX_INDEX(fend_x.levels[mi], fend_y.levels[mi], msizes[mi])];
	} else {
	  cerr << "Error: unknown scope" << endl;
	  exit(1);
	}

        total += value;
      }
  return total;
}

int main(int argc, char **argv)
{
  UserParams params;
  parse_user_arguments(argc, argv, params);
  cerr << params;

  // list of model features
  vector < vector<double> > model_matrices(params.m_features.size);

  int dim_x = 1;
  for (int i=0; i<params.o_features_x.size; i++)
    dim_x *= params.o_features_x.to[i] - params.o_features_x.from[i] + 1;
  int dim_y = 1;
  for (int i=0; i<params.o_features_y.size; i++)
    dim_y *= params.o_features_y.to[i] - params.o_features_y.from[i] + 1;

  // container of output bins
  map<string, Bin> obins_x;
  map<string, Bin> obins_y;

  // read fends file
  cerr << "Reading fends into memory ..." << endl;
  read_fends_table(params.fends_fn, params.m_features, params.o_features_x, params.o_features_y,
		   params.from_fend, params.to_fend, obins_x, obins_y);

  // read model feature files
  if (params.model != m_count)
    read_feature_tables(params.m_features.fns, params.m_features.sizes, model_matrices);

  // open output file
  ofstream out(params.ofn.c_str());
  massert(out.is_open(), "could not open file %s", params.ofn.c_str());

  out << join_strings(params.o_features_x.fields.begin(), params.o_features_x.fields.end(), "\t", NumberToString(1)) << "\t";
  out << join_strings(params.o_features_y.fields.begin(), params.o_features_y.fields.end(), "\t", NumberToString(2)) << "\t";
  out << "value\n";

  int total_bins = dim_x * dim_y;
  cerr << "Traversing all fend pairs, number of bins: " << total_bins << endl;
  int count = 0;

  vector<string> keys_x;
  vector<string> keys_y;
  init_keys(params.o_features_x, 0, string(""), keys_x);
  init_keys(params.o_features_y, 0, string(""), keys_y);

  vector<int>& msizes = params.m_features.sizes;

  for (unsigned int ix=0; ix<keys_x.size(); ix++)
    for (unsigned int iy=0; iy<keys_y.size(); iy++)
      {
        map<string, Bin>::iterator it_x = obins_x.find(keys_x[ix]);
        map<string, Bin>::iterator it_y = obins_y.find(keys_y[iy]);

        if (it_x == obins_x.end() || it_y == obins_y.end()) {
	  //	  out << keys_x[ix] << "\t" << keys_y[iy] << "\t" << 0 << endl;
	  continue;
	}

        vector<Fend>& fends_x = (*it_x).second.fends;
        vector<Fend>& fends_y = (*it_y).second.fends;

        vector<Fend> fends_x_unique;
        vector<Fend> fends_y_unique;
        vector<Fend> fends_common;
        intersect_fends(fends_x, fends_y, fends_x_unique, fends_y_unique, fends_common);

        double total_common = fends_expected(fends_common, fends_common, params.scope, params.model,
					     model_matrices, msizes,
					     params.prior);
        double total = fends_expected(fends_x, fends_y, params.scope, params.model,
				      model_matrices, msizes,
				      params.prior);
	total -= total_common/2;

	if (total > 0) {
	  char buffer[100];
	  sprintf(buffer, "%g", total);
	  out << keys_x[ix] << "\t" << keys_y[iy] << "\t" << buffer << endl;
	}

        // progress report
        if (++count % 100 == 0)
	  fprintf(stderr, "%d/%d (%.1f%%)\n", count, total_bins, 100.0 * count / total_bins);
      }
  fprintf(stderr, "%d/%d (%.1f%%)\n", count, total_bins, 100.0 * count / total_bins);

  out.close();
  cerr << "done" << endl;

  return 0;
}
