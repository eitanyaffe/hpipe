#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <cmath>
#include <time.h>

#include <queue>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////
// Utils
//////////////////////////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures
//////////////////////////////////////////////////////////////////////////////////////////////////

enum ScoreType { stContacts, stShared, stSharedCorrected, stSharedCorrectedReg };
/*
  Clustering metric:
   stContacts: log2 of number of contacts between contigs plus 1 (log2(C_12+1)
   stShared: number of shared neighbours of the two contigs (S_12)
   stSharedCorrected: shared neighbours corrected for expected shared function of marignals:
       (1/N + S_12) / ((M1*M2)/N)
   where N is total contigs and M1/2 is contig marginals
   stSharedCorrectedReg: same as stSharedCorrected but regulated differently to prefer contigs with many share
       (1 + S_12) / (1 + (M1*M2)/N)
*/

ScoreType str2ScoreType(string str) {
  if (str == "contacts")
    return stContacts;
  if (str == "shared")
    return stShared;
  if (str == "shared_corrected")
    return stSharedCorrected;
  if (str == "shared_corrected_reg")
    return stSharedCorrectedReg;
  cout << "unknown score type: " << str << endl;
  exit(-1);
}

struct Contig
{
  // string iden
  string name;

  // length in bp
  int length;

  // index
  int index;

  // containing node
  int node_index;

  // set of contig neighbours
  set<int> contig_neighbors;

  // number of contacts
  int n_neighbors;

  // clustering order
  int order;
};

// node in cluster tree
struct Node
{
  // index
  int index;

  // all contained contigs
  vector<Contig*> contigs;

  // node children
  bool leaf;
  int left_i, right_i, parent_i;

  // nodes were united based on this score
  double score;

  // flag when node was already united
  bool root_flag;

  // map from node to contact count
  map<int, double> node_neighbors;

  // clustering order
  int order;

  // visual params
  int level;
  int start_contig, end_contig;
};

struct ScorePrint
{
  string contig1, contig2;
  double score;
  ScorePrint(string _c1, string _c2, double _s) :
    contig1(_c1), contig2(_c2), score(_s) {};
};

// read contigs from input file
vector<Contig> g_contigs;

// contact matrix, from contig indices to count
map<int, map<int, int> > g_matrix;

Contig& get_contig(unsigned int index) {
  massert(index >= 0 && index < g_contigs.size(), "contig index out of bounds");
  return (g_contigs[index]);
}

double get_score(Contig& c1, Contig& c2, ScoreType stype)
{
  // first treat stContacts, which are different all together
  if (stype == stContacts)
    return log2(g_matrix[c1.index][c2.index]+1);

  const int N = g_contigs.size();

  vector<int> v;
  set<int>& s1 = c1.contig_neighbors;
  set<int>& s2 = c2.contig_neighbors;
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(v , v.end()));

  double obs = v.size();
  double exp = double(c1.n_neighbors * c2.n_neighbors) / double(N);
  massert(exp > 0, "expected must be > 0");

  double result = 0;
  switch(stype) {
  case stShared:
    result = obs;
    break;
  case stSharedCorrected:
    result = log2((obs+1.0/double(N)) / exp);
    break;
  case stSharedCorrectedReg:
    result = log2((obs+1) / (exp+1));
    break;
  default:
    cout << "unknown score type" << endl;
    exit(-1);
  }

  return (result);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Parsing and handling user arguments
//////////////////////////////////////////////////////////////////////////////////////////////////

class UserParams {
public:
  string contig_ifn; // contig table
  string contact_ifn; // contact table
  int min_coverage;
  int ratio;
  string ofn;
  ScoreType stype;

  // ctor with default values
  UserParams() {
    ofn = "score.out";
    contig_ifn = "i/c";
    contact_ifn = "i/m";
    min_coverage = 2;
    stype = stContacts;
  }
};

void usage(const char* name, UserParams& params)
{
  fprintf(stderr, "usage: %s [options]\n", name);
  cout << " -contigs <fn>: contig table" << endl;
  cout << " -contacts <fn>: contact table" << endl;
  cout << " -min_coverage <int>: minimal number of adjacent contigs per contig" << endl;
  cout << " -ratio <int>: sample all contig pairs with this ratio" << endl;
  cout << " -score_type <str>: clustering metric" << endl;
  cout << " -out <fn>: output" << endl;
  cout << "Clustering metric: " << endl;
  cout << "   contacts: log of number of contacts between contigs + 1 (log2(C_12+1)" << endl;
  cout << "   shared: number of shared neighbours of the two contigs (S_12)" << endl;
  cout << "   shared_corrected: shared neighbours corrected for expected shared function of marignals: " << endl;
  cout << "       (1/N + S_12) / ((M1*M2)/N)" << endl;
  cout << "   where N is total contigs and M1/2 is contig marginals" << endl;
  cout << "   shared_corrected_reg: same as stSharedCorrected but regulated differently to prefer contigs with many share" << endl;
  cout << "       (1 + S_12) / (1 + (M1*M2)/N)" << endl;

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

      if (option == "-contigs")
	params.contig_ifn = arg;
      else if (option == "-contacts")
	params.contact_ifn = arg;
      else if (option == "-ratio")
	params.ratio = atoi(arg);
      else if (option == "-min_coverage")
	params.min_coverage = atof(arg);
      else if (option == "-score_type")
	params.stype = str2ScoreType(arg);
      else if (option == "-out")
	params.ofn = arg;
      else {
	cout << "Error: unknown option: " << option << endl;
	exit(1);
      }

      i += 2;
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

//////////////////////////////////////////////////////////////////////////////////////////////////
// I/O utility functions
//////////////////////////////////////////////////////////////////////////////////////////////////

void split_line(istream &in, vector<string> &fields, char delim)
{
  fields.resize(0);
  string field;
  while(in) {
    char c = in.get();
    if (c == -1)
      break;
    if(c == '\r') {
      continue;
    }
    if(c == '\n') {
      fields.push_back(field);
      field.resize(0);
      break;
    }
    if(c == delim) {
      fields.push_back(field);
      field.resize(0);
    } else {
      field.push_back(c);
    }
  }

  if (field.length() > 0)
    fields.push_back(field);
}

int get_field_index(string field, const vector<string>& titles)
{
  int result = -1;
  for (unsigned int i = 0; i<titles.size(); i++)
    if (titles[i] == field)
      result = i;
  if (result == -1) {
    cout << "unknown field " << field << endl;
    exit(1);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// I/O
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_contigs(string fn, vector<Contig>& contigs, map<string, int>& contig_map, map<string, set<string> >& marginals, int min_coverage)
{
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  int index = 0;
  // parse header
  split_line(in, fields, delim);
  int id_idx = get_field_index("contig", fields);
  int length_idx = get_field_index("length", fields);
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string name = fields[id_idx];

    if ((int)marginals[name].size() < min_coverage)
      continue;

    Contig contig;
    contig.name = name;
    contig.length = (length_idx != -1) ? atoi(fields[length_idx].c_str()) : 0;
    contig.index = index++;
    contig.node_index = -1;
    contig.n_neighbors = 0;
    contig.order = -1;
    contigs.push_back(contig);

    // save map from name to index for contact table
    contig_map[contig.name] = contig.index;
  }
  in.close();
}

int read_contacts(string fn, vector<pair< string, string> >& contacts, map<string, set<string> >& marginals)
{
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  int result = 0;

  // parse header
  split_line(in, fields, delim);
  int contig1_idx = get_field_index("contig1", fields);
  int contig2_idx = get_field_index("contig2", fields);
  int contacts_idx = get_field_index("contacts", fields);
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string contig1 = fields[contig1_idx];
    string contig2 = fields[contig2_idx];

    // skip self contacts
    if (contig1 == contig2)
      continue;

    int n_contacts = atoi(fields[contacts_idx].c_str());
    if (n_contacts < 1)
      continue;

    contacts.push_back(make_pair(contig1, contig2));
    marginals[contig1].insert(contig2);
    marginals[contig2].insert(contig1);
    result++;
  }
  in.close();

  return (result);
}

int process_contacts(map<string, int>& contig_map, vector<Contig>& contigs,
		     vector<pair< string, string> >& contacts,
		     vector<ScorePrint>& scores, ScoreType stype,
		     map<int, map<int, int> >& matrix)
{
  cout << "going over contig pairs ..." << endl;
  int result = 0;
  for (unsigned int i=0; i < contacts.size(); i++) {

    string contig1_str = contacts[i].first;
    string contig2_str = contacts[i].second;

    // skip contigs if not found in table
    if ( (contig_map.find(contig1_str) == contig_map.end()) ||
	 (contig_map.find(contig2_str) == contig_map.end()) )
      continue;

    result++;

    int contig_index1 = contig_map[contig1_str];
    int contig_index2 = contig_map[contig2_str];

    Contig& contig1 = contigs[contig_index1];
    Contig& contig2 = contigs[contig_index2];

    contig1.contig_neighbors.insert(contig_index2);
    contig2.contig_neighbors.insert(contig_index1);

    // update matrix
    matrix[contig_index1][contig_index2]++;
    matrix[contig_index2][contig_index1]++;
  }

  // connect contig to itself by definition
  for( unsigned int i=0; i < contigs.size(); i++)
    contigs[i].contig_neighbors.insert(contigs[i].index);

  for( unsigned int i=0; i < contigs.size(); i++)
    contigs[i].n_neighbors = contigs[i].contig_neighbors.size();

  return (result);
}

void output_results(UserParams& params, ScoreType stype, int ratio)
{
  map<double, int> result;
  int N = g_contigs.size();
  cout << "total number of pairs: " << N*N << endl;
  cout << "sample ratio: " << ratio << endl;
  for (int i1 = 0; i1 < N; i1++) {
  for (int i2 = 0; i2 < N; i2++) {
    if (i1 % ratio != 0)
      continue;
    Contig& c1 = get_contig(i1);
    Contig& c2 = get_contig(i2);
    double score = get_score(c1, c2, stype);
    result[score]++;
  } }

  cout << "result has " << result.size() << " bins" << endl;

  cout << "saving metric table: " << params.ofn << endl;
  // output contigs
  ofstream out(params.ofn.c_str());
  out << "score\tcount" << endl;
  massert(out.is_open(), "could not open file %s", params.ofn.c_str());
  for (map<double, int>::iterator it = result.begin(); it != result.end(); ++it)
    out << (*it).first << "\t" << (*it).second << endl;
  out.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  UserParams params;
  parse_user_arguments(argc, argv, params);

  // first read contacts
  cout << "reading contact table: " << params.contact_ifn << endl;
  vector<pair< string, string> > contacts;
  map<string, set<string> > marginals;
  int n_contacts = read_contacts(params.contact_ifn, contacts, marginals);
  massert(n_contacts > 0, "no contact pairs found");
  cout << "total number of contigs: " << marginals.size() << endl;
  cout << "number of contacts: " << n_contacts << endl;

  // map from contig string identifier to contig index
  map<string, int> contig_map;

  // init contigs and contacts
  cout << "contig table: " << params.contig_ifn << endl;
  read_contigs(params.contig_ifn, g_contigs, contig_map, marginals, params.min_coverage);
  massert(g_contigs.size() > 0, "no contigs pass filter");
  cout << "number of contigs with a min marginal coverage >" << params.min_coverage << " is: " << g_contigs.size() << endl;

  // set contig and node connectivity
  vector<ScorePrint> scores;
  process_contacts(contig_map, g_contigs, contacts, scores, params.stype, g_matrix);

  // output metric
  output_results(params, params.stype, params.ratio);

  return (0);
}
