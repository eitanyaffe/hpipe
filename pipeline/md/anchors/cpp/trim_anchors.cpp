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

#include "Model.h"
#include "util.h"
#include "Params.h"

using namespace std;

struct ContigAnchorPair
{
  int m_contig_index;
  int m_anchor_index;
  double m_score;
  int m_version;
  ContigAnchorPair(int contig_index, int anchor_index, double score, int version) :
    m_contig_index(contig_index), m_anchor_index(anchor_index), m_score(score), m_version(version) {};
};

class CompareContigAnchorPair
{
public:
  CompareContigAnchorPair() {};
  bool operator() (const ContigAnchorPair& lhs, const ContigAnchorPair& rhs) const
  {
    return (lhs.m_score < rhs.m_score);
  }
};

typedef priority_queue<ContigAnchorPair, vector<ContigAnchorPair>, CompareContigAnchorPair> Queue;

enum TrimType { ttMulti, ttLowScore, ttFewObserved, ttSmallAnchor, ttNormal };

string TrimType2string(TrimType ttype) {
  switch (ttype) {
  case ttMulti: return string("multi");
  case ttLowScore: return string("lowScore");
  case ttFewObserved: return string("fewObserved");
  case ttSmallAnchor: return string("smallAnchor");
  case ttNormal: return string("normal");
  default: mexit("unknown type: %d", (int)ttype);
  }
  return NULL;
}

class AnchorTrimmer
{
private:
  vector<double> m_threshold_exp;
  vector<int> m_threshold_obs;

  int m_min_contacts;
  int m_min_anchor_size;
  Model m_model;

  // number of contigs
  int m_num_contigs;

  // number of ancnors
  int m_num_anchors;

  // removed contigs
  int m_num_contigs_trimmed;

  map<int, map<int, int> > m_contacts;
  vector<int> m_original_anchor;
  vector<int> m_offensive_anchor;
  vector<TrimType> m_trimtype;

  map<int, string> m_fend_contig_map;
  map<string, int> m_contig_index_map;
  map<int, string> m_rev_contig_index_map;

  vector< vector<int> > m_obs_matrix;
  vector< vector<double> > m_exp_matrix;

  vector<Contig> m_contigs;
  vector< set <int > > m_anchors;
  vector< int > m_anchor_size;

  vector<int> m_contig_length;

  map<pair<int, int>, int > m_version_map;
  Queue m_queue;

  void remove_contig(int contig_index, int contig_anchor, TrimType ttype, int violating_anchor=-1, double score=0);
  void update_pair(int contig_index, int anchor_index, int obs, double exp);

  void inspect_anchor(int anchor);


public:
  AnchorTrimmer(int _min_contacts, int _min_anchor_size, double prior, ModelFeatures features) :
    m_min_contacts(_min_contacts), m_min_anchor_size(_min_anchor_size), m_model(features, prior) {};

  inline bool is_significant_contact(double exp, int obs);

  // io
  void read_contacts(string fn);
  void read_ca_matrix(string fn);
  void read_model(string fn);
  void read_threshold_table(string fn);
  void read_contig_length(string fn);

  void init();
  void verify_counts();
  void trim();

  void write_ca_matrix(string fn);
  void save(string ofn, string ofn_matrix);
};

void init_params(int argc, char **argv, Parameters& params)
{

  // filenames
  params.add_parser("fends", new ParserFilename("input binned fend table"), true);
  params.add_parser("contacts", new ParserFilename("input fend observed pairs"), true);
  params.add_parser("contigs", new ParserFilename("table with contig/length columns"), true);
  params.add_parser("ca_matrix", new ParserFilename("contig-anchor start matrix"), true);
  params.add_parser("threshold_table", new ParserFilename("file with threshold table"), true);
  params.add_parser("prior_fn", new ParserFilename("prior fn"), true);
  params.add_parser("ofn", new ParserFilename("output file", "o"), false);
  params.add_parser("ofn_matrix", new ParserFilename("output matrix file", "om"), false);

  // optional parameters
  params.add_parser("verify_ca_matrix", new ParserBoolean("should verify contig-anchor start matrix", false), false);
  params.add_parser("min_contacts", new ParserInteger("discard multi with at least min_contacts between contig and anchor", 10), false);
  params.add_parser("min_anchor_size", new ParserInteger("min_anchor_size", 10), false);

  // model
  params.add_parser("model_num", new ParserInteger("number of model parameters"), true);
  params.add_parser("f_model_field_i", new ParserString("model field name i", "", true), false);
  params.add_parser("f_model_fn_i", new ParserFilename("model file i", "", true), false);
  params.add_parser("f_model_size_i", new ParserInteger("model matrix dimension i", 0, true), false);

  if (argc == 1) {
    params.usage(argv[0]);
    exit(0);
  }

  // read command line params
  params.read(argc, argv);
  params.parse(true);

  int N = params.get_int("model_num");
  for (int i=1; i<=N; i++) {
    string s = to_string((long long)i);
    params.add_parser("f_model_size_" + s, new ParserInteger("model matrix dimension "+ s), true);
    params.add_parser("f_model_fn_" + s, new ParserFilename("model file " + s), true);
    params.add_parser("f_model_field_" + s, new ParserString("model field name " + s), true);
  }
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

void AnchorTrimmer::read_model(string fn)
{
  m_model.read_fends_table(fn, "fend", "contig", "coord", "anchor", m_num_anchors);
  cout << "number of anchors: " << m_num_anchors << endl;

  m_contig_index_map = m_model.contig_index_map();
  m_rev_contig_index_map = m_model.rev_contig_index_map();
  m_fend_contig_map = m_model.fend_contig_map();

  m_contigs = m_model.contigs();
  m_num_contigs = m_contigs.size();
  cout << "number of contigs: " << m_num_contigs << endl;
}

void AnchorTrimmer::init()
{
  m_anchors.resize(m_num_anchors);
  m_anchor_size.resize(m_num_anchors);

  m_original_anchor.resize(m_num_contigs);
  m_offensive_anchor.resize(m_num_contigs);
  m_trimtype.resize(m_num_contigs);

  m_num_contigs_trimmed = 0;

  for (int i=0; i<m_num_anchors; i++)
    m_anchor_size[i] = 0;

  for (int i=0; i<m_num_contigs; i++) {
    massert(m_contigs[i].anchor >= 0 && m_contigs[i].anchor < m_num_anchors, "invalid anchor: %d", m_contigs[i].anchor);
    m_anchors[m_contigs[i].anchor].insert(i);
    m_original_anchor[i] = m_contigs[i].anchor;
    m_offensive_anchor[i] = -1;
    m_trimtype[i] = ttNormal;

    // compute anchor size
    m_anchor_size[m_contigs[i].anchor] += m_contig_length[i];
  }
}

void AnchorTrimmer::read_contacts(string fn)
{
  cout << "reading contacts file: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int fend1_idx = get_field_index("fend1", fields);
  int fend2_idx = get_field_index("fend2", fields);

  int total = 0;
  int atotal = 0;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    total++;
    int fend1 = atoi(fields[fend1_idx].c_str());
    int fend2 = atoi(fields[fend2_idx].c_str());
    if (m_fend_contig_map.find(fend1) == m_fend_contig_map.end() || m_fend_contig_map.find(fend2) == m_fend_contig_map.end())
      continue;

    if ((m_contig_index_map.find(m_fend_contig_map[fend1]) == m_contig_index_map.end()) ||
	(m_contig_index_map.find(m_fend_contig_map[fend2]) == m_contig_index_map.end()))
      continue;

    int contig1 = m_contig_index_map[m_fend_contig_map[fend1]];
    int contig2 = m_contig_index_map[m_fend_contig_map[fend2]];

    // skip self contacts
    if (contig1 == contig2)
      continue;

    atotal++;
    m_contacts[contig1][contig2] += 1;
    m_contacts[contig2][contig1] += 1;
  }
  in.close();
  cout << "total contacts: " << total << endl;
  cout << "anchor contacts: " << atotal << endl;
}

void AnchorTrimmer::read_contig_length(string fn)
{
  m_contig_length.resize(m_contig_index_map.size());

  cout << "reading contigs table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int contig_idx = get_field_index("contig", fields);
  int length_idx = get_field_index("length", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string contig_id = fields[contig_idx].c_str();
    int length = atoi(fields[length_idx].c_str());

    if (m_contig_index_map.find(contig_id) != m_contig_index_map.end()) {
      int contig_index = m_contig_index_map[contig_id];
      m_contig_length[contig_index] = length;
    }
  }
  in.close();
}

void AnchorTrimmer::read_threshold_table(string fn)
{
  cout << "reading threshold table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int exp_idx = get_field_index("exp", fields);
  int obs_idx = get_field_index("obs", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    int obs = atoi(fields[obs_idx].c_str());
    double exp = atof(fields[exp_idx].c_str());

    m_threshold_exp.push_back(exp);
    m_threshold_obs.push_back(obs);
  }
}

void AnchorTrimmer::write_ca_matrix(string fn)
{
 cout << "writing ca matrix file: " << fn << endl;
  ofstream out(fn.c_str());
  massert(out.is_open(), "could not open file %s", fn.c_str());

  out << "contig\tanchor\toriginal_anchor\tobserved\texpected" << endl;
  for (int j=0; j<m_num_anchors; j++) {
    vector<int>& anchor_obs = m_obs_matrix[j];
    vector<double>& anchor_exp = m_exp_matrix[j];
    for (int i=0; i<m_num_contigs; i++) {
      int obs = anchor_obs[i];
      double exp = anchor_exp[i];

      if ((obs == 0) || (exp == 0))
	continue;
      string contig = m_rev_contig_index_map[i];
      out << contig << "\t" << j+1 << "\t" << m_original_anchor[i]+1 << "\t" << obs << "\t" << exp << endl;
    }
  }
  out.close();
}

void AnchorTrimmer::read_ca_matrix(string fn)
{
  cout << "reading ca matrix file: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int contig_idx = get_field_index("contig", fields);
  int anchor_idx = get_field_index("other_cell", fields);
  int obs_idx = get_field_index("observed_contact_count", fields);
  int exp_idx = get_field_index("expected", fields);

  m_obs_matrix.resize(m_num_anchors);
  m_exp_matrix.resize(m_num_anchors);
  for (int j=0; j<m_num_anchors; j++) {
    m_obs_matrix[j].resize(m_num_contigs);
    m_exp_matrix[j].resize(m_num_contigs);
  }

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string contig_id = fields[contig_idx];
    if (m_contig_index_map.find(contig_id) == m_contig_index_map.end())
      continue;

    int contig = m_contig_index_map[contig_id];
    massert(contig >= 0 && contig < m_num_contigs, "contig %d out of range %d", contig, m_num_contigs);

    int anchor = atoi(fields[anchor_idx].c_str()) - 1;
    int obs = atoi(fields[obs_idx].c_str());
    double exp = atof(fields[exp_idx].c_str());

    massert(anchor >= 0 && anchor < m_num_anchors, "anchor %d out of range %d", anchor, m_num_anchors);

    m_obs_matrix[anchor][contig] = obs;
    m_exp_matrix[anchor][contig] = exp;

  }
  in.close();
}

void AnchorTrimmer::verify_counts()
{
  cout << "verfiy anchor-contig matrices" << endl;
  for (int j=0; j<m_num_anchors; j++) {
    vector<int>& anchor_obs = m_obs_matrix[j];
    vector<double>& anchor_exp = m_exp_matrix[j];
    set<int>& anchor_contigs = m_anchors[j];

    cout << "contig count of anchor " << j+1 << ": " << anchor_contigs.size() << endl;

    for (int i=0; i<m_num_contigs; i++) {
      int obs = anchor_obs[i];
      int v_obs = 0;
      for (set<int>::iterator ik = anchor_contigs.begin(); ik != anchor_contigs.end(); ++ik) {
	if (i != *ik)
	  v_obs += m_contacts[i].find(*ik) != m_contacts[i].end() ? m_contacts[i][*ik] : 0;
      }

      if (v_obs != obs) {
	printf("contig %s, anchor %d\n", m_rev_contig_index_map[i].c_str(), j+1);
	printf("pre-computed obs: %d != contacts file obs: %d\n", obs, v_obs);
	exit(-1);
      }

      double exp = anchor_exp[i];
      double v_exp = m_model.compute(m_contigs[i].id, j);
      if (abs(v_exp-exp)/exp > 0.001) {
	printf("contig %s, anchor %d\n", m_rev_contig_index_map[i].c_str(), j+1);
	printf("pre-computed exp: %f != model exp: %f\n", exp, v_exp);
	exit(-1);
      }
    }
  }
}

void AnchorTrimmer::remove_contig(int contig_index, int contig_anchor, TrimType ttype, int violating_anchor, double score)
{
  if (contig_anchor == -1)
    return;

  cout << "removed contig: "
       << "N=" << m_num_contigs_trimmed+1
       << ", type=" << TrimType2string(ttype)
       << ", contig=" << m_rev_contig_index_map[contig_index].c_str()
       << ", anchor=" << contig_anchor+1
       << ", anchor_bp=" << m_anchor_size[contig_anchor]
       << ", secondary anchor=" << violating_anchor+1
       << ", score=" << score << endl;

  if (ttype != ttMulti) {
    massert(violating_anchor == -1, "if not ttMulti v anchor must be -1");
  } else {
    massert(violating_anchor != -1, "if ttMulti v anchor must not be -1");
  }
  m_anchors[contig_anchor].erase(contig_index);
  m_contigs[contig_index].anchor = -1;
  m_offensive_anchor[contig_index] = violating_anchor;
  m_trimtype[contig_index] = ttype;
  m_num_contigs_trimmed++;

  m_anchor_size[contig_anchor] -= m_contig_length[contig_index];

  // update for all other contigs the contig-anchor obs and exp matrices and if needed add to queue
  vector<int>& anchor_obs = m_obs_matrix[contig_anchor];
  vector<double>& anchor_exp = m_exp_matrix[contig_anchor];
  for (int i=0; i<m_num_contigs; i++) {
    if (i == contig_index)
      continue;
    anchor_obs[i] -= (m_contacts[i].find(contig_index) != m_contacts[i].end()) ? m_contacts[i][contig_index] : 0;
    anchor_exp[i] -= m_model.compute(m_contigs[i].id, m_contigs[contig_index].id);

    if (anchor_obs[i] < 0) {
      cout << "Error: obs below zero, contig_index=" << contig_index << ", i=" << i
	   << ", o=" << m_contacts[i][contig_index]
	   << ", contig_index=" << m_rev_contig_index_map[contig_index].c_str()
	   << ", violating_anchor=" << violating_anchor+1
	   << ", contig_i=" << m_rev_contig_index_map[i].c_str()
	   << ", anchor_i=" << m_contigs[i].anchor+1 << endl;
      exit(-1);
    }
    if (anchor_exp[i] < 0)
      anchor_exp[i] = 0;
  }
}

void AnchorTrimmer::inspect_anchor(int anchor)
{
  if (m_anchors[anchor].size() == 0)
    return;

  vector<int>& anchor_obs = m_obs_matrix[anchor];
  vector<double>& anchor_exp = m_exp_matrix[anchor];

  while (1) {
    set<int>& anchor_contigs = m_anchors[anchor];
    int min_contig = -1;

    // trim: ttSmallAnchor
    if (m_anchor_size[anchor] < m_min_anchor_size) {
      set<int> tmp = m_anchors[anchor];
      for (set<int>::iterator it = tmp.begin(); it != tmp.end(); ++it)
	remove_contig((*it), anchor, ttSmallAnchor);
      break;
    }

    // trim: ttFewObserved
    int min_obs = 0;
    for (set<int>::iterator it = anchor_contigs.begin(); it != anchor_contigs.end(); ++it) {
      int contig_index = (*it);
      int obs = anchor_obs[contig_index];
      if (obs < m_min_contacts && ((min_contig == -1) || (obs < min_obs))) {
	min_contig = contig_index;
	min_obs = obs;
      }
    }
    if (min_contig != -1) {
      remove_contig(min_contig, anchor, ttFewObserved);
      continue;
    }

    // trim: ttLowScore
    double local_min_score = 0;
    for (set<int>::iterator it = anchor_contigs.begin(); it != anchor_contigs.end(); ++it) {
      int contig_index = (*it);
      int obs = anchor_obs[contig_index];
      double exp = anchor_exp[contig_index];
      if (exp == 0) {
	min_contig = contig_index;
	break;
      }

      double score = log(obs/exp);
      if (((min_contig == -1) || (score < local_min_score)) && !is_significant_contact(exp, obs)) {
	min_contig = contig_index;
	local_min_score = score;
      }
    }
    if (min_contig != -1)
      remove_contig(min_contig, anchor, ttLowScore);

    if (min_contig == -1)
      break;
  }
}

bool AnchorTrimmer::is_significant_contact(double exp, int obs)
{
  int low_index = lower_bound (m_threshold_exp.begin(), m_threshold_exp.end(), exp) - m_threshold_exp.begin() - 1;
  if (low_index >= (int)m_threshold_exp.size()) low_index = (int)m_threshold_exp.size() - 1;
  if (low_index < 0) low_index = 0;
  int obs_t = m_threshold_obs[low_index];
  return (obs > obs_t);

}

void AnchorTrimmer::update_pair(int contig_index, int anchor_index, int obs, double exp)
{
  int contig_anchor = m_contigs[contig_index].anchor;
  string contig = m_rev_contig_index_map[contig_index];

  // !!!
  //  string contig = m_rev_contig_index_map[contig_index];
  // if (contig == string("c3734"))
  //  cout << "anchor=" << anchor_index << ", contig=" << contig << endl;

  if (m_anchors[anchor_index].size() == 0)
    return;

  pair<int,int> key = make_pair(anchor_index, contig_index);
  int version = (m_version_map.find(key) != m_version_map.end()) ? m_version_map[key] + 1 : 0;
  bool added = false;
  if ((exp > 0) && (contig_anchor != anchor_index)) {
    // schedule for removal inter with a high score
    double score = log(obs/exp);

    //    cerr << "contig=" << contig << ", contig_anchor=" << contig_anchor << ", anchor=" << anchor_index << ", obs=" << obs << ", exp=" << exp << ", score=" << score << ", is_significant=" << (is_significant_contact(exp, obs) ? "T" : "F") << endl;

    if (is_significant_contact(exp, obs)) {
      m_queue.push(ContigAnchorPair(contig_index, anchor_index, score, version));
      added = true;
    }
  }

  // invalidate old pair if a new one added or the queue already contained a prior version
  if (version>0 || added)
    m_version_map[make_pair(anchor_index,contig_index)] = version;
}

void AnchorTrimmer::save(string ofn, string ofn_matrix)
{
  cout << "number of removed contigs: " << m_num_contigs_trimmed << endl;

  cout << "writing file: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  out << "contig\tanchor\tdiscarded\ttype\tsecondary_anchor" << endl;

  for (unsigned int i=0; i<m_contigs.size(); i++)
    out << m_contigs[i].id << "\t" << m_original_anchor[i]+1 << "\t"
	<< (m_contigs[i].anchor == -1 ? "T" : "F") << "\t"
	<< TrimType2string(m_trimtype[i]) << "\t"
	<< m_offensive_anchor[i]+1 << endl;
  out.close();

  // write final matrix
  write_ca_matrix(ofn_matrix);
}

void AnchorTrimmer::trim()
{
  // first clear all weak connections to original anchors and small anchors
  cout << "pre-inspecting all anchors..." << endl;
  for (int i=0; i<m_num_anchors; i++)
    inspect_anchor(i);
  cout << "number of contig-anchor pairs with few supporting contacts trimmed: " << m_num_contigs_trimmed << endl;

  // add relevant anchor-contig pairs to queue
  cout << "initializing queue" << endl;
  for (int j=0; j<m_num_anchors; j++) {
    vector<int>& anchor_obs = m_obs_matrix[j];
    vector<double>& anchor_exp = m_exp_matrix[j];
    for (int i=0; i<m_num_contigs; i++)
      update_pair(i, j, anchor_obs[i], anchor_exp[i]);
  }

  cout << "starting queue size: " << m_queue.size() << endl;
  while (m_queue.size() > 0) {
    // find top violating contig with max o/e which is above threshold
    ContigAnchorPair top = m_queue.top();
    m_queue.pop();

    int removed_contig = top.m_contig_index;
    int violating_anchor = top.m_anchor_index;
    int original_anchor = m_contigs[removed_contig].anchor;
    double score = top.m_score;
    int version = top.m_version;

    if (m_contigs[removed_contig].anchor == -1)
      continue;

    massert(original_anchor != violating_anchor, "violating and original anchors must be different");

    // skip invalidated ContigAnchorPair
    pair<int,int> key = make_pair(violating_anchor, removed_contig);
    massert(m_version_map.find(key) != m_version_map.end(), "contig-anchor pair not found in version map of size: %d", m_version_map.size());
    if (version < m_version_map[key])
      continue;

    // remove contig from anchor
    remove_contig(removed_contig, original_anchor, ttMulti, violating_anchor, score);
    inspect_anchor(original_anchor);

    vector<int>& anchor_obs = m_obs_matrix[original_anchor];
    vector<double>& anchor_exp = m_exp_matrix[original_anchor];
    for (int i=0; i<m_num_contigs; i++) {
      if (i == removed_contig)
	continue;
      update_pair(i, original_anchor, anchor_obs[i], anchor_exp[i]);
    }
  }

  if ((int)m_contigs.size() == m_num_contigs_trimmed) {
    cout << "no contigs left, adjust trimming parameters" << endl;
    exit(-1);
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  Parameters params;
  init_params(argc, argv, params);

  ModelFeatures features; // model features
  features.size = params.get_int("model_num");
  for (int i=1; i<=features.size; i++) {
    string s = to_string((long long)i);
    features.fields.push_back(params.get_string("f_model_field_" + s));
    features.fns.push_back(params.get_string("f_model_fn_" + s));
    features.sizes.push_back(params.get_int("f_model_size_" + s));
  }

  int min_contacts = params.get_int("min_contacts");
  string table_fn = params.get_string("threshold_table");
  string fends_fn = params.get_string("fends");
  string contigs_fn = params.get_string("contigs");
  string contacts_fn = params.get_string("contacts");
  string ca_matrix_fn = params.get_string("ca_matrix");
  bool verify_ca_matrix = params.get_bool("verify_ca_matrix");
  int min_anchor_size = params.get_int("min_anchor_size");
  string ofn = params.get_string("ofn");
  string ofn_matrix = params.get_string("ofn_matrix");

  double prior = read_double(params.get_string("prior_fn"));
  cout << "prior: " << prior << endl;


  AnchorTrimmer atrimmer(min_contacts, min_anchor_size, prior, features);

  // read threshold table
  atrimmer.read_threshold_table(table_fn);

  // model
  atrimmer.read_model(fends_fn);

  // contacts[i][j] := number of contacts between contig i and contig j
  atrimmer.read_contacts(contacts_fn);

  // read contig length
  atrimmer.read_contig_length(contigs_fn);

  atrimmer.init();

  cout << "initializing anchor-contig matrices" << endl;
  atrimmer.read_ca_matrix(ca_matrix_fn);

  if (verify_ca_matrix)
    atrimmer.verify_counts();

  atrimmer.trim();

  atrimmer.save(ofn, ofn_matrix);
}
