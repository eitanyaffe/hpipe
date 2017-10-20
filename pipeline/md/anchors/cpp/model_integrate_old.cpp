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

  2. A list of 0 or more model features. Each model feature is defined by:
  - field name: same name as the one in the fend table
  - file: contains correction matrix for feature
  - size: matrix dimension

  3. Two list (for dimension x and y) of 1 or more output features. Each feature is defined by:
  - field name: same name as the one in the fend table
  - from: feature start index
  - to: feature to index

  4. A prior on the interaction probablity

  Output
  ------
  The program outputs a single file that contains expected counts (according to model)
  for the requested output bins.

  Examples
  --------
  The following counts the size of the bins:
  %> expected_count fends.table 1 0 1 coord_bin 1 1000 1 coord_bin 1 1000

  The following checks model performance (fragment length only) in gc plane:
  %> expected_count fends.table 1 1 frag_len_bin frag_len.f 50 1 gc_bin 1 50 1 gc_bin 1 50

  Remarks
  -------
  All input and output indices are 1-based. The model matrix is represented as a vector,
  so it is zero based.

*/
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures
//////////////////////////////////////////////////////////////////////////////////////////////////


struct CisTables {
  // one decay file per chrom
  vector< vector<double> > cis_table_uniform;

  // for two way clustering use 5 decay profiles
  vector< vector<double> > cis_table_intra_1;
  vector< vector<double> > cis_table_intra_2;
  vector< vector<double> > cis_table_inter_1;
  vector< vector<double> > cis_table_inter_2;
  vector< vector<double> > cis_table_mixed;
};

// The model has a set of features. Each feature is factor which is a function of
// the two interacting fends
struct ModelFeatures
{
  int size;
  vector<string> fields;
  vector<string> fns;
  vector<int> sizes;
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
  int chrom_code;
  bool is_chrom_na;
  double original_coord, virtual_coord;
  int id;
  vector<int> levels;
  int cluster;
  int domain;

  string inter_id;
};

// used for compaction purposes
struct Frag
{
  int index;
  int chrom_code;
  int coord;
  double frag_length;
  double lfactor;
};

// every bin holds a collection of all fends that have the same levels for all model features
struct Bin
{
  vector<Fend> fends;
};

// filter on contact type
enum Filter { f_both, f_trans, f_close_cis, f_far_cis, f_id, f_inter };
Filter string2filter(string str)
{
  Filter filter;
  if (str == "both")
    filter = f_both;
  else if (str == "trans")
    filter = f_trans;
  else if (str == "close_cis")
    filter = f_close_cis;
  else if (str == "far_cis")
    filter = f_far_cis;
  else if (str == "id")
    filter = f_id;
  else if (str == "inter")
    filter = f_inter;
  else {
    cerr << "Unknown filter " << str << endl;
    exit(1);
  }
  return (filter);
}

string filter2string(Filter filter)
{
  string result;
  switch (filter) {
  case f_both:
    result = "both";
    break;
  case f_trans:
    result = "trans";
    break;
  case f_close_cis:
    result = "close_cis";
    break;
  case f_far_cis:
    result = "far_cis";
    break;
  case f_id:
    result = "id";
    break;
  case f_inter:
    result = "inter";
    break;
  default:
    cerr << "Unknown filter index " << filter << endl;
    exit(1);
  }
  return (result);
}

enum CisMode { cm_none, cm_uniform, cm_2clusters};
CisMode string2cismode(string str)
{
  CisMode cismode;
  if (str == "none")
    cismode = cm_none;
  else if (str == "uniform")
    cismode = cm_uniform;
  else if (str == "2clusters")
    cismode = cm_2clusters;
  else {
    cerr << "Unknown cismode " << str << endl;
    exit(1);
  }
  return (cismode);
}

string cismode2string(CisMode cismode)
{
  string result;
  switch (cismode) {
  case cm_none:
    result = "none";
    break;
  case cm_uniform:
    result = "uniform";
    break;
  case cm_2clusters:
    result = "2clusters";
    break;
  default:
    cerr << "Unknown cismode, index=" << cismode << endl;
    exit(1);
  }
  return (result);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Parsing and handling user arguments
//////////////////////////////////////////////////////////////////////////////////////////////////

class UserParams {
public:
  string fends_fn; // input file (.binned)
  string ofn; // output file (.e_contact)

  int from_fend, to_fend; // limit output to fend range
  Filter filter; // trans/cis filter on output
  int cis_threshold; // threshold between close and far cis

  // cis parameters, relevant only if normalizing with cis
  CisMode cis_mode; // none / uniform / 2clusters
  double cis_log_step; // cis table binsize (log10 scale)
  string cis_table_directory; // dir of cis tables

  // relevant for compaction
  string compaction_fn; // if specified we correct for compaction as well

  double trans_prior; // if not normalizing with cis we use the trans prior for contact

  ModelFeatures m_features; // model features

  // inter id field, relevant only for filter=inter
  string inter_field;

  // output features
  OutputFeatures o_features_x;
  OutputFeatures o_features_y;

  // ctor with default values
  UserParams() {
    fends_fn = "results/wtb_s0_trans.binned";
    ofn = "out";
    from_fend = 0;
    to_fend = 0;
    filter = f_both;
    cis_threshold = -1;
    cis_mode = cm_none;
    cis_table_directory = "";
    cis_log_step = 0;
    compaction_fn = "";
    trans_prior = 1;
    inter_field = "";
  }
};

ostream& operator <<(ostream &os, const UserParams &params)
{
  os << "binned fends table: " << params.fends_fn << endl;
  os << "output file: " << params.ofn << endl;

  os << "fend from: " << params.from_fend << endl;
  os << "fend to: " << params.to_fend << endl;

  os << "cismode: " << cismode2string(params.cis_mode) << endl;
  if (params.cis_mode == cm_none) {
    os << "not using cis correction, only trans prior is used for all contacts" << endl;
  } else {
    os << "cis table log step: " << params.cis_log_step << endl;
    os << "cis table directory: " << params.cis_table_directory << endl;
    if (params.cis_mode == cm_uniform) {
      os << "normalizing with uniform cis decay" << endl;
    } else if (params.cis_mode == cm_2clusters) {
      os << "normalizing with 2-cluster cis decay" << endl;
    } else {
      cerr << "Error: unknown cis mode" << endl;
      exit(1);
    }
  }

  os << "filter: " << filter2string(params.filter) << endl;
  os << "cis threshold: " << params.cis_threshold << endl;

  if (params.filter == f_inter)
    os << "filter inter id field: " << params.inter_field << endl;

  if (params.compaction_fn != "")
    os << "compaction ifn: " << params.compaction_fn << endl;
  else
    os << "not using compaction parameters" << endl;

  os << "trans_prior: " << params.trans_prior << endl;

  // model features
  if (params.m_features.size == 0)
    os << "No model functions" << endl;
  else for (int i=0; i<params.m_features.size; i++)
	 os << "M" << i << ": " << params.m_features.fields[i] << ", " <<
	   params.m_features.fns[i] << ", dim=" << params.m_features.sizes[i] << endl;

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
  cerr << " -cis_mode <none|uniform|2clusters>: how to handle cis" << endl;
  cerr << " -cis_log_step <log10(bp)>: cis table binsize" << endl;
  cerr << " -cis_table_dir <fn>: directory of all cis tables" << endl;
  cerr << " -compaction_fn <fn>: input compaction table" << endl;
  cerr << " -from_fend <fend index>: fend mask start" << endl;
  cerr << " -to_fend <fend index>: fend mask end" << endl;
  cerr << " -filter <type>: output only both|trans|close_cis|far_cis|id|inter contacts" << endl;
  cerr << " -inter_field <string>: field defining inter fend relations, relevant for inter filter" << endl;
  cerr << " -cis_threshold <bp>: threshold between close and far cis" << endl;
  cerr << " -model_num <n>: number of model features" << endl;
  cerr << " -f_model_field <string>: model feature field" << endl;
  cerr << " -f_model_fn <filename>: model feature filename" << endl;
  cerr << " -f_model_size <n>: model feature matrix dimension" << endl;
  cerr << " -output_{x|y}_num <n>: number of output x|y features" << endl;
  cerr << " -f_{x|y}_field <string>: output x|y feature field" << endl;
  cerr << " -f_{x|y}_from <int>: output x|y from bin index" << endl;
  cerr << " -f_{x|y}_to <int>: output x|y to bin index" << endl;
  cerr << " -output <fn>: output file" << endl;

  cerr << "Example:\n" <<
    "%> ./bin/model_integrate -fends results/wtb_s0_trans_10k.cbinned -cis_log_step 0.1 -cis_table_dir results/cis_decays/ -compaction_fn results/wtb_s0_trans_50_200.cmp -from_fend 0 -to_fend 0 -filter both -cis_threshold 100000 -model_num 2 -f_model_field frag_gc_bin -f_model_fn results/wtb_s0_trans_frag_gc_bin.f -f_model_size 20 -f_model_field frag_len_bin -f_model_fn results/wtb_s0_trans_frag_len_bin.f -f_model_size 20 -output_x_num 1 -f_x_field coord_bin -f_x_from 1 -f_x_to 100 -output_y_num 1 -f_y_field coord_bin -f_y_from 1 -f_y_to 100 -output out.e_contact" << endl;
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
      else if (option == "-inter_field")
	params.inter_field = string(arg);
      else if (option == "-from_fend")
	params.from_fend = atoi(arg);
      else if (option == "-to_fend")
	params.to_fend = atoi(arg);
      else if (option == "-cis_mode")
	params.cis_mode = string2cismode(string(arg));
      else if (option == "-cis_table_dir")
	params.cis_table_directory = arg;
      else if (option == "-cis_log_step")
	params.cis_log_step = atof(arg);
      else if (option == "-filter")
	params.filter = string2filter(string(arg));
      else if (option == "-cis_threshold")
	params.cis_threshold = atoi(arg);
      else if (option == "-compaction_fn")
	params.compaction_fn = arg;
      else if (option == "-trans_prior")
	params.trans_prior = atof(arg);
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
  if (params.m_features.size != (int)params.m_features.fns.size()) {
    cerr << "Error: Expecting " << params.m_features.size << " model feature filenames" << endl;
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

// use a global structure for chromosome codes
static map<string, int> chrom_codes;
static map<int, string> rev_chrom_codes;
static int next_code = 0;
int chr2int(string chr)
{
  if (chrom_codes.find(chr) == chrom_codes.end())
    {
      chrom_codes[chr] = next_code;
      rev_chrom_codes[next_code] = chr;
      next_code++;
    }
  return chrom_codes[chr];
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
    exit(-1);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read fragment compaction file
//////////////////////////////////////////////////////////////////////////////////////////////////

struct comp_frags : public binary_function<Frag, Frag, bool> {
  bool operator()(Frag x, Frag y) {
    if (x.chrom_code != y.chrom_code)
      return (x.chrom_code < y.chrom_code);
    else
      return (x.coord < y.coord);
  }
};

void read_compaction_table(const string& fn, map <int, pair<double, double> >& frag_coords_map)
{
  // map from chromosome code to vector of fragments
  map<int, vector<Frag> > chr_frags;

  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int frag_ind = get_field_index("frag", fields, fn);
  int chr_ind = get_field_index("chr", fields, fn);
  int coord_ind = get_field_index("coord", fields, fn);
  int flen_ind = get_field_index("len", fields, fn);
  int lfactor_ind = get_field_index("lfactor", fields, fn);

  // so we don't add a frag twice
  map<int, bool> frag_added_map;

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    int index = atoi(fields[frag_ind].c_str());
    int chrom_code = chr2int(fields[chr_ind]);

    if (frag_added_map.find(index) != frag_added_map.end())
      continue;
    frag_added_map[index] = true;

    Frag frag;
    frag.index = index;
    frag.coord = atoi(fields[coord_ind].c_str());
    frag.frag_length = atof(fields[flen_ind].c_str());
    frag.lfactor = atof(fields[lfactor_ind].c_str());

    // add fragment
    chr_frags[chrom_code].push_back(frag);
  }
  in.close();

  // go over chroms, compute virtual coords for all fragments
  map<int, vector<Frag> >::iterator it = chr_frags.begin();
  while (it != chr_frags.end())
    {
      vector<Frag>& frags = (*it).second;
      // sort fragments for each chromosome
      sort(frags.begin(), frags.end(), comp_frags());
      //         double coord = 0;
      //         for (int i = 0; i < (int)frags.size(); i++)
      //         {
      //             Frag& frag = frags[i];
      //             double flen = frag.frag_length;
      //             pair<double, double> coords(coord, coord+flen);
      //             frag_coords_map[frag.index] = coords;
      //             coord += flen;
      //         }

      double total_offset = 0;
      for (int i = 0; i < (int)frags.size(); i++)
        {
	  Frag& frag = frags[i];
	  double flen = frag.frag_length;
	  double lfactor = frag.lfactor;
	  // offset = flen - original_flen = flen - flen / factor = flen - flen/(2^lfactor))
	  pair<double, double> coords(frag.coord+total_offset, frag.coord+total_offset+flen);
	  frag_coords_map[frag.index] = coords;

	  // update cumulative offset
	  double offset = flen * (1-1/pow(2,lfactor));
	  total_offset += offset;

	  //printf("frag=%d coord=%d len=%f lfactor=%f offset=%f total_offset=%f\n",
	  //       frag.index, frag.coord, flen, lfactor, offset, total_offset);
        }
      it++;
    }

  // debug: dump fragment map
  //     cout << "frag\tstart\tend" << endl;
  //     for (map <int, pair<double, double> >::iterator it = frag_coords_map.begin();
  //              it != frag_coords_map.end(); it++)
  //         cout << (*it).first << "\t" << (*it).second.first << "\t" << (*it).second.second << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read fend file
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_fends_table(const string& fn, const ModelFeatures& m_features,
                      const OutputFeatures& o_features_x, const OutputFeatures& o_features_y,
                      const int from_fend, const int to_fend,
                      map <int, pair<double, double> >& frag_coords_map,
                      map<string, Bin>& obins_x, map<string, Bin>& obins_y, Filter filter,
                      bool use_clusters, bool use_domains, string inter_field)
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
  int frag_ind = get_field_index("frag", fields, fn);
  int chr_ind = get_field_index("contig", fields, fn);
  int strand_ind = get_field_index("strand", fields, fn);
  int coord_ind = get_field_index("coord", fields, fn);
  int cluster_ind = use_clusters ? get_field_index("cluster", fields, fn) : 0;
  int domain_ind = use_domains ? get_field_index("domain", fields, fn) : 0;

  int id_ind = -1;
  if (filter == f_id)
    id_ind = get_field_index("id", fields, fn);

  int inter_ind = -1;
  if (filter == f_inter)
    inter_ind = get_field_index(inter_field, fields, fn);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      return;

    Fend fend;
    fend.index = atoi(fields[fend_ind].c_str());
    fend.chrom_code = chr2int(fields[chr_ind]);
    fend.is_chrom_na = (fields[chr_ind] == "0");
    fend.cluster = use_clusters ? atoi(fields[cluster_ind].c_str()) : 1;
    if (fend.cluster != 1 && fend.cluster != 2) {
      cerr << "Error: currently supports two clusters only" << endl;
      exit(1);
    }

    // use domain if required
    fend.domain = use_domains ? atoi(fields[domain_ind].c_str()) : -1;

    // keep original coord for cis threshold
    fend.original_coord = atof(fields[coord_ind].c_str());

    // coord can depend on fragment compactions if compaction was supplied
    if (frag_coords_map.size() > 0) {
      int frag_index = atoi(fields[frag_ind].c_str());
      if (frag_coords_map.find(frag_index) == frag_coords_map.end()) {
	cerr << "Error: fragment " << frag_index << " not defined" << endl;
	exit(1);
      }
      fend.virtual_coord = (fields[strand_ind] == "+") ?
	frag_coords_map[frag_index].first :
	frag_coords_map[frag_index].second;
    } else {
      fend.virtual_coord = fend.original_coord;
    }

    fend.id = -1;
    if (filter == f_id)
      fend.id = atoi(fields[id_ind].c_str());
    if (filter == f_trans && fend.is_chrom_na)
      continue;

    fend.inter_id = -1;
    if (filter == f_inter)
      fend.inter_id = fields[inter_ind];

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
// Read cis table
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_cis_table(string fn, double binsize, double& table_log_start, vector<double>& result)
{
  cerr << "reading cis table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';
  bool first = true;

  split_line(in, fields, delim);
  int bin_ind = get_field_index("log_dist", fields, fn);
  int prob_ind = get_field_index("ratio", fields, fn);
  double min_prob = 1;
  double prev_log_dist = 0;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    double log_dist = atof(fields[bin_ind].c_str());
    double prob = atof(fields[prob_ind].c_str());

    // keep first logdist, if there are multiple tables all must start with the same bin
    if (first && table_log_start != 0 && table_log_start != log_dist) {
      cerr << "cis tables have different starting bins: " << table_log_start << ", " << log_dist << endl;
      exit(1);
    }
    first = false;
    if (table_log_start == 0)
      table_log_start = log_dist;

    // verify binsize
    if (prev_log_dist != 0) {
      double file_binsize = log_dist - prev_log_dist;
      if (fabs(file_binsize - binsize) > 0.001) {
	fprintf(stderr, "binsize in table (%f) does not match user-specified binsize (%f)\n",
		file_binsize, binsize);
	exit(1);
      }
    }

    result.push_back(prob);
    if (min_prob > prob && prob > 0)
      min_prob = prob;
    prev_log_dist = log_dist;
  }
  in.close();

  double tail_count = 2;
  for (unsigned int i = 0; i < result.size(); i++) {
    if (result[i] == 0)
      result[i] = min_prob / tail_count++;
  }
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
	  cerr << " " << fend.index << " " << fend.chrom_code << " " << rev_chrom_codes[fend.chrom_code] << ": ";
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
                      Filter filter, const double cis_threshold,
                      const vector < vector<double> >& model_matrices, const vector<int>& msizes,
                      const double trans_prior, const CisMode cis_mode,
                      const double cis_table_log_step, const double cis_table_log_start, const double cis_table_min_dist,
                      const CisTables& cis_tables)
{
  double total = 0;
  vector<Fend>::const_iterator jt_x, jt_y;

  for(jt_x = fends_x.begin(); jt_x != fends_x.end(); jt_x++)
    for(jt_y = fends_y.begin(); jt_y != fends_y.end(); jt_y++)
      {
        const Fend& fend_x = *jt_x;
        const Fend& fend_y = *jt_y;
        //cout << fend_x.chrom_code << " " << fend_y.chrom_code << endl;
        switch ( filter ) {
        case f_id:
	  if (fend_x.id != fend_y.id)
	    continue;
	  break;
        case f_inter:
	  if (fend_x.inter_id == fend_y.inter_id)
	    continue;
	  break;
        case f_trans:
	  if ((fend_x.chrom_code == fend_y.chrom_code) || fend_x.is_chrom_na || fend_y.is_chrom_na)
	    continue;
	  break;
        case f_close_cis:
	  if (fend_x.chrom_code != fend_y.chrom_code || (fabs(fend_x.original_coord - fend_y.original_coord) > cis_threshold))
	    continue;
	  break;
        case f_far_cis:
	  if (fend_x.chrom_code != fend_y.chrom_code || (fabs(fend_x.original_coord - fend_y.original_coord) < cis_threshold))
	    continue;
	  break;
        case f_both:
	  break;
        default:
	  cerr << "Error: unknown filter" << endl;
	  exit(1);
	}

        // skip if self-interaction
        if (fend_x.index == fend_y.index)
	  continue;

        double value = trans_prior;

        if ((fend_x.chrom_code == fend_y.chrom_code) && (cis_mode != cm_none)) {
	  const vector<double>* cis_table = NULL;
	  int chrom_code = fend_x.chrom_code;

	  // determine which decay table to use
	  if (cis_mode == cm_uniform)
	    cis_table = &cis_tables.cis_table_uniform[chrom_code];
	  else if (cis_mode == cm_2clusters) {
	    assert(fend_x.domain != -1 && fend_y.domain != -1);
	    if (fend_x.domain == fend_y.domain) {
	      // intra
	      assert(fend_x.cluster == fend_y.cluster);
	      if (fend_x.cluster == 1)
		cis_table = &cis_tables.cis_table_intra_1[chrom_code];
	      else
		cis_table = &cis_tables.cis_table_intra_2[chrom_code];
	    } else {
	      // inter
	      if (fend_x.cluster == fend_y.cluster) {
		if (fend_x.cluster == 1)
		  cis_table = &cis_tables.cis_table_inter_1[chrom_code];
		else
		  cis_table = &cis_tables.cis_table_inter_2[chrom_code];
	      } else {
		cis_table = &cis_tables.cis_table_mixed[chrom_code];
	      }
	    }
	  }

	  if (cis_table == NULL || cis_table->size() == 0) {
	    cerr << "Error: cis table not defined" << endl;
	    exit(1);
	  }

	  // here we use virtual coord which might take compaction into account
	  double coord_dist = fabs(fend_x.virtual_coord - fend_y.virtual_coord);
	  double dist = (coord_dist < cis_table_min_dist) ? 0 : (log10(coord_dist) - cis_table_log_start) / cis_table_log_step;
	  double before_dist = floor(dist);
	  double after_dist = ceil(dist);

	  // correct if exactly integer
	  if (after_dist == before_dist)
	    after_dist += 1;

	  double w_before = after_dist - dist;
	  double w_after = dist - before_dist;
	  int bin_before = (int)before_dist;
	  int bin_after = (int)after_dist;

	  int cis_table_last_index = (int)cis_table->size()-1;
	  assert(bin_after >= 0);
	  if (bin_before > cis_table_last_index)
	    bin_before = cis_table_last_index;
	  if (bin_after > cis_table_last_index)
	    bin_after = cis_table_last_index;

	  double before_value = (bin_before >= 0) ? (*cis_table)[bin_before] : 1;
	  value = w_before*before_value + w_after * (*cis_table)[bin_after];
        }

        for(unsigned int mi = 0; mi < model_matrices.size(); mi++)
	  value *= model_matrices[mi][MATRIX_INDEX(fend_x.levels[mi], fend_y.levels[mi], msizes[mi])];

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

  // if compaction_fn is specified we compute virtual coords for all fragments
  map <int, pair<double, double> > frag_coords_map; // map from fragment to a pair of corrected coordinates (start,end)
  if (params.compaction_fn != "") {
    cerr << "Reading compaction table into memory ..." << endl;
    read_compaction_table(params.compaction_fn, frag_coords_map);
  }

  // read fends file
  bool use_clusters = (params.cis_mode == cm_2clusters);
  bool use_domains = (params.cis_mode == cm_2clusters);
  cerr << "Reading fends into memory ..." << endl;
  read_fends_table(params.fends_fn, params.m_features, params.o_features_x, params.o_features_y,
		   params.from_fend, params.to_fend, frag_coords_map, obins_x, obins_y, params.filter, use_clusters, use_domains, params.inter_field);

  // read cis tables if needed
  CisTables cis_tables;
  double cis_table_log_start = 0;
  int n_chroms = rev_chrom_codes.size();
  if (params.cis_mode == cm_uniform) {
    cis_tables.cis_table_uniform.resize(n_chroms);
    for (int i = 0; i < n_chroms; i++) {
      assert(rev_chrom_codes.find(i) != rev_chrom_codes.end());
      string chr = rev_chrom_codes[i];
      string cis_table_fn = params.cis_table_directory + "/" + chr + ".cis_decay";
      read_cis_table(cis_table_fn, params.cis_log_step, cis_table_log_start, cis_tables.cis_table_uniform[i]);
    }
  }
  if (params.cis_mode == cm_2clusters) {
    cis_tables.cis_table_intra_1.resize(n_chroms);
    cis_tables.cis_table_intra_2.resize(n_chroms);
    cis_tables.cis_table_inter_1.resize(n_chroms);
    cis_tables.cis_table_inter_2.resize(n_chroms);
    cis_tables.cis_table_mixed.resize(n_chroms);
    for (int i = 0; i < n_chroms; i++) {
      assert(rev_chrom_codes.find(i) != rev_chrom_codes.end());
      string chr = rev_chrom_codes[i];

      // intra decay
      string cis_table_1_intra_fn = params.cis_table_directory + "/" + chr + "_0_intra.cis_decay";
      string cis_table_2_intra_fn = params.cis_table_directory + "/" + chr + "_1_intra.cis_decay";

      // use strict inter decay
      string cis_table_1_inter_fn = params.cis_table_directory + "/" + chr + "_0_inter_strict.cis_decay";
      string cis_table_2_inter_fn = params.cis_table_directory + "/" + chr + "_1_inter_strict.cis_decay";

      // for inter combos use the mixed
      string cis_table_mixed_fn = params.cis_table_directory + "/" + chr + "_mixed.cis_decay";

      read_cis_table(cis_table_1_intra_fn, params.cis_log_step, cis_table_log_start, cis_tables.cis_table_intra_1[i]);
      read_cis_table(cis_table_2_intra_fn, params.cis_log_step, cis_table_log_start, cis_tables.cis_table_intra_2[i]);
      read_cis_table(cis_table_1_inter_fn, params.cis_log_step, cis_table_log_start, cis_tables.cis_table_inter_1[i]);
      read_cis_table(cis_table_2_inter_fn, params.cis_log_step, cis_table_log_start, cis_tables.cis_table_inter_2[i]);
      read_cis_table(cis_table_mixed_fn, params.cis_log_step, cis_table_log_start, cis_tables.cis_table_mixed[i]);
    }
  }
  double cis_table_min_dist = pow(10, cis_table_log_start);

  // read model feature files
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

        if (it_x == obins_x.end() || it_y == obins_y.end())
	  {
            out << keys_x[ix] << "\t" << keys_y[iy] << "\t" << 0 << endl;
            continue;
	  }

        vector<Fend>& fends_x = (*it_x).second.fends;
        vector<Fend>& fends_y = (*it_y).second.fends;

        vector<Fend> fends_x_unique;
        vector<Fend> fends_y_unique;
        vector<Fend> fends_common;
        intersect_fends(fends_x, fends_y, fends_x_unique, fends_y_unique, fends_common);

        double total_common = fends_expected(fends_common, fends_common, params.filter, params.cis_threshold,
					     model_matrices, msizes, params.trans_prior,
					     params.cis_mode, params.cis_log_step, cis_table_log_start, cis_table_min_dist, cis_tables);
        double total = fends_expected(fends_x, fends_y, params.filter, params.cis_threshold,
				      model_matrices, msizes, params.trans_prior,
				      params.cis_mode, params.cis_log_step, cis_table_log_start, cis_table_min_dist, cis_tables);
	total -= total_common/2;

        char buffer[100];
        sprintf(buffer, "%g", total);
        if (total > 0)
	  out << keys_x[ix] << "\t" << keys_y[iy] << "\t" << buffer << endl;

        // progress report
        if (++count % 100 == 0)
	  fprintf(stderr, "%d/%d (%.1f%%)\n", count, total_bins, 100.0 * count / total_bins);
      }
  fprintf(stderr, "%d/%d (%.1f%%)\n", count, total_bins, 100.0 * count / total_bins);

  out.close();
  cerr << "done" << endl;

  return 0;
}
