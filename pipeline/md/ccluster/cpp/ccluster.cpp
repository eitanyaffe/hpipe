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

class CompareNode
{
public:
  CompareNode() {};
  bool operator() (const Node* lhs, const Node* rhs) const
  {
    return (lhs->contigs.size() < rhs->contigs.size());
  }
};

typedef priority_queue<Node*, vector<Node*>, CompareNode> NodeQueue;

struct NodePair
{
  int node_index1, node_index2;
  double score;
  NodePair(int i1, int i2, double s) : node_index1(i1), node_index2(i2), score(s) {};
};

class CompareNodePair
{
public:
  CompareNodePair() {};
  bool operator() (const NodePair& lhs, const NodePair& rhs) const
  {
    return (lhs.score < rhs.score);
  }
};

typedef priority_queue<NodePair, vector<NodePair>, CompareNodePair> NodePairQueue;

// read contigs from input file
vector<Contig> g_contigs;

// dynamic list of cluster nodes
vector<Node> g_nodes;

// cache scores
map < pair<int, int>, double> g_contig_score_cache;

// contact matrix, from contig indices to count
map<int, map<int, int> > g_matrix;

Node& get_node(unsigned int index)
{
  massert(index >= 0 && index < g_nodes.size(), "node index out of bounds");
  return (g_nodes[index]);
}

Node& add_node(double score)
{
  int new_index = g_nodes.size();
  g_nodes.resize(g_nodes.size() + 1);
  Node& node = get_node(new_index);
  node.index = new_index;

  node.score = score;
  node.leaf = false;
  node.root_flag = true;

  return (node);
}

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

double compute_edge_score(Node& n1, Node& n2, map<int, double>& node_neighbors_1, int max_contig_pairs, ScoreType stype)
{
  if (node_neighbors_1.find(n2.index) != node_neighbors_1.end())
    return (node_neighbors_1[n2.index]);

  int n1_count = n1.contigs.size();
  int n2_count = n2.contigs.size();

  int score_count = 0;
  double score_sum = 0;
  if (max_contig_pairs == -1 || (n1_count*n2_count < max_contig_pairs)) {
    score_count = n1.contigs.size()* n2.contigs.size();
    for (unsigned int i1 = 0; i1 < n1.contigs.size(); i1++)
      for (unsigned int i2 = 0; i2 < n2.contigs.size(); i2++)
	score_sum += get_score(*n1.contigs[i1], *n2.contigs[i2], stype);
  } else {
    score_count = max_contig_pairs;
    for (int i=0; i<max_contig_pairs; i++) {
      unsigned int i1 = rand() % n1_count;
      unsigned int i2 = rand() % n2_count;
      score_sum += get_score(*n1.contigs[i1], *n2.contigs[i2], stype);
    }
  }

  return score_sum / (double)score_count;
}

int update_neighbors(Node& p, Node& l, Node& r, Node& n, double min_score, int max_contig_pairs, ScoreType stype, bool correct_cluster_size)
{
  // cout << "p=" << p.index << ", l=" << l.index << ", r=" << r.index << ", n=" << n.index << endl;

  // weights childrens
  double N_l = l.contigs.size();
  double N_r = r.contigs.size();
  double N = N_l + N_r;

  // score of children
  double l_score = compute_edge_score(l, n, l.node_neighbors, max_contig_pairs, stype);
  double r_score = compute_edge_score(r, n, r.node_neighbors, max_contig_pairs, stype);

  double w_l = N_l/N;
  double w_r = N_r/N;
  double score = w_l * l_score + w_r * r_score;

  // add to neighbor list only if score is high enough
  if ((min_score == -1) || score/(correct_cluster_size ? N : 1) > min_score) {
    p.node_neighbors[n.index] = score;
    n.node_neighbors[p.index] = score;
    return 1;
  } {
    return 0;
  }
}

int merge_nodes(Node& p, Node& l, Node& r, double min_score, int max_contig_pairs, ScoreType stype, bool correct_cluster_size)
{
  //  cout << "merging nodes: (" << l.index << "," << r.index << ")" << endl;
  massert(l.index != r.index, "cannot unite node with itself");

  p.left_i = l.index;
  p.right_i = r.index;

  l.root_flag = false;
  r.root_flag = false;

  l.parent_i = p.index;
  r.parent_i = p.index;

  // concat contig vectors which don't overlap
  p.contigs = l.contigs;
  p.contigs.insert(p.contigs.end(), r.contigs.begin(), r.contigs.end());

  // number of new node edges
  int count = 0;

  // handle neighbors of left child
  map<int, bool> visited;
  for (map<int, double>::iterator it = l.node_neighbors.begin(); it != l.node_neighbors.end(); ++it) {
    if (visited[it->first] || (it->first == l.index) || (it->first == r.index))
      continue;
    count += update_neighbors(p, l, r, get_node(it->first), min_score, max_contig_pairs, stype, correct_cluster_size);
    visited[it->first] = true;
  }

  // handle neighbors of right child
  for (map<int, double>::iterator it = r.node_neighbors.begin(); it != r.node_neighbors.end(); ++it) {
    if (visited[it->first] || (it->first == l.index) || (it->first == r.index))
      continue;
    count += update_neighbors(p, l, r, get_node(it->first), min_score, max_contig_pairs, stype, correct_cluster_size);
    visited[it->first] = true;
  }

  return count;
}

void init_queue(NodePairQueue& queue, double min_score, bool correct_cluster_size)
{
  int pcount = 0;
  for (vector<Node>::iterator it = g_nodes.begin(); it != g_nodes.end(); ++it) {
    Node& node1 = *it;
    for (map<int,double>::iterator jt=node1.node_neighbors.begin();
	 jt != node1.node_neighbors.end(); ++jt) {
      Node& node2 = get_node(jt->first);
      if (node1.index >= node2.index)
	continue;

      double factor = correct_cluster_size ? (node1.contigs.size() + node2.contigs.size()) : 1;
      double score = jt->second;

      if ((min_score != -1) && (score/factor < min_score)) continue;

      NodePair node_pair(node1.index, node2.index, score);
      queue.push(node_pair);
      pcount++;
    }
  }
  // cout << "initial pair candidates: " << pcount << endl;
}

int update_queue(NodePairQueue& queue, Node& new_node, double min_score, bool correct_cluster_size)
{
  int pcount = 0;
  for (map<int,double>::iterator it=new_node.node_neighbors.begin(); it != new_node.node_neighbors.end(); ++it) {
    Node& node_o = get_node(it->first);

    if (!node_o.root_flag || (new_node.index == node_o.index)) continue;

    double factor = correct_cluster_size ? (node_o.contigs.size() + new_node.contigs.size()) : 1;
    double score = it->second;
    if ((min_score != -1) && (score/factor < min_score)) continue;

    pcount++;
    NodePair node_pair(new_node.index, node_o.index, score);
    queue.push(node_pair);
  }
  return pcount;
  //  cout << "new pair candidates: " << pcount << endl;
}

void init_nodes(vector<Contig>& contigs, vector<Node>& nodes)
{
  // we reserve enough space for a full binary tree
  nodes.reserve(contigs.size() * 2 - 1);

  // start off by one node per contig
  nodes.resize(contigs.size());

  for (unsigned int i=0; i<contigs.size(); i++) {
    Contig* contig = &contigs[i];
    Node& node = nodes[i];
    node.index = i;
    node.leaf = true;
    node.root_flag = true;
    node.score = -1;
    node.contigs.resize(0);
    node.contigs.push_back(contig);
    contig->node_index = node.index;
  }
}

void dfs_tree(Node& node, vector<Node>& nodes, int& contig_order, int& node_order, int level)
{
  node.order = node_order++;
  node.level = level;
  if (node.leaf) {
    massert(node.contigs.size() == 1, "should have exactly 1 contig in leaf node");
    Contig* contig = node.contigs[0];
    contig->order = contig_order++;
    node.start_contig = contig->order;
    node.end_contig = contig->order+1;
  } else {
    Node& left = nodes[node.left_i];
    Node& right = nodes[node.right_i];

    dfs_tree(left, nodes, contig_order, node_order, level+1);
    dfs_tree(right, nodes,  contig_order, node_order, level+1);

    node.start_contig = min(left.start_contig, right.start_contig);
    node.end_contig = max(left.end_contig, right.end_contig);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Parsing and handling user arguments
//////////////////////////////////////////////////////////////////////////////////////////////////

class UserParams {
public:
  string contig_ifn; // contig table
  string contact_ifn; // contact table

  string contig_ofn;
  string tree_ofn;

  string score_ofn;

  int min_coverage;
  double min_score;

  bool correct_cluster_size;

  int max_contig_pairs;

  ScoreType stype;

  // ctor with default values
  UserParams() {
    score_ofn = "score.out";
    contig_ifn = "i/c";
    contact_ifn = "i/m";
    max_contig_pairs = -1;
    min_score = 0;
    min_coverage = 2;
    correct_cluster_size = false;
    contig_ofn = "contigs.out";
    tree_ofn = "tree.out";
    stype = stContacts;
  }
};

void usage(const char* name, UserParams& params)
{
  fprintf(stderr, "usage: %s [options]\n", name);
  cout << " -contigs <fn>: contig table" << endl;
  cout << " -contacts <fn>: contact table" << endl;
  cout << " -max_contig_pairs <int>: sample score when number of pairs is larger than this threshold" << endl;
  cout << " -min_score <double>: stop when reaching below this score" << endl;
  cout << " -cluster_normalize <T|F>: min_score corrected for cluster size" << endl;
  cout << " -min_coverage <int>: minimal number of adjacent contigs per contig" << endl;
  cout << " -score_type <str>: clustering metric" << endl;
  cout << " -out_contigs <fn>: contig table output" << endl;
  cout << " -out_tree <fn>: tree output" << endl;
  cout << " -out_scores <fn>: score output" << endl;
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
      else if (option == "-max_contig_pairs")
	params.max_contig_pairs = atoi(arg);
      else if (option == "-min_score")
	params.min_score = atof(arg)/100;
      else if (option == "-cluster_normalize")
	params.correct_cluster_size = string(arg) == "T";
      else if (option == "-min_coverage")
	params.min_coverage = atof(arg);
      else if (option == "-score_type")
	params.stype = str2ScoreType(arg);
      else if (option == "-out_contigs")
	params.contig_ofn = arg;
      else if (option == "-out_scores")
	params.score_ofn = arg;
      else if (option == "-out_tree")
	params.tree_ofn = arg;
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
		     vector<Node>& nodes, vector<pair< string, string> >& contacts,
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

  cout << "creating nodes..." << endl;
  for (unsigned int i=0; i < contacts.size(); i++) {

    string contig1_str = contacts[i].first;
    string contig2_str = contacts[i].second;

    // skip contigs if not found in table
    if ( (contig_map.find(contig1_str) == contig_map.end()) ||
	 (contig_map.find(contig2_str) == contig_map.end()) )
      continue;

    int contig_index1 = contig_map[contig1_str];
    int contig_index2 = contig_map[contig2_str];

    Contig& contig1 = contigs[contig_index1];
    Contig& contig2 = contigs[contig_index2];

    // finally connect nodes
    int node_index1 = contig1.node_index;
    int node_index2 = contig2.node_index;

    massert( (node_index1 != -1) && (node_index2 != -1), "nodes not defined");

    double score = get_score(contig1, contig2, stype);
    nodes[node_index1].node_neighbors[node_index2] = score;
    nodes[node_index2].node_neighbors[node_index1] = score;

    scores.push_back(ScorePrint(contig1_str, contig2_str, score));
  }

  return (result);
}

void print_node(ofstream& out, Node& node, vector<Node>& nodes)
{
  double left_score = (node.leaf) ? -1 : nodes[node.left_i].score;
  double right_score = (node.leaf) ? -1 : nodes[node.right_i].score;
  double parent_score = (node.root_flag) ? -1 : nodes[node.parent_i].score;
  out << node.index << "\t" << node.order << "\t" << node.level << "\t" << node.leaf << "\t"
      << node.contigs.size() << "\t"
      << node.start_contig << "\t"
      << node.end_contig << "\t"
      << node.score << "\t"
      << node.parent_i << "\t" << node.left_i << "\t" << node.right_i << "\t"
      << parent_score << "\t" << left_score << "\t" << right_score << endl;

  if (!node.leaf) {
    Node& left = nodes[node.left_i];
    Node& right = nodes[node.right_i];
    print_node(out, left, nodes);
    print_node(out, right, nodes);
  }
}

void output_results(UserParams& params, vector<int>& roots)
{
  // set contig order
  int contig_order = 1;
  int node_order = 1;
  for (unsigned int i=0; i<roots.size(); i++) {
    Node& node = get_node(roots[i]);
    dfs_tree(node, g_nodes, contig_order, node_order, 0);
  }

  cout << "saving contig table: " << params.contig_ofn << endl;

  // output contigs
  ofstream out(params.contig_ofn.c_str());
  massert(out.is_open(), "could not open file %s", params.contig_ofn.c_str());
  int N = g_contigs.size();
  out << "contig\tlength\torder\tnode\tmarginal\n";
  for (int i = 0; i < N; i++) {
    Contig& contig = get_contig(i);
    out << contig.name << "\t" << contig.length << "\t" << contig.order
	<< "\t" << contig.node_index << "\t" << contig.n_neighbors
	<< endl;
  }
  out.close();

  cout << "saving tree table: " << params.tree_ofn << endl;

  // output tree
  ofstream out_tree(params.tree_ofn.c_str());
  massert(out_tree.is_open(), "could not open file %s", params.tree_ofn.c_str());
  out_tree << "node\torder\tlevel\tleaf\ttree_size\tstart_contig\tend_contig\tscore\tparent\tleft\tright\tparent_score\tleft_score\tright_score\n";
  for (unsigned int i=0; i<roots.size(); i++) {
    Node& node = get_node(roots[i]);
    print_node(out_tree, node, g_nodes);
  }
  out_tree.close();
}

void output_scores(UserParams& params, vector<ScorePrint>& scores)
{
  cout << "saving score table: " << params.score_ofn << endl;

  // output contigs
  ofstream score_out(params.score_ofn.c_str());
  massert(score_out.is_open(), "could not open file %s", params.score_ofn.c_str());
  score_out << "contig1\tcontig2\tscore\n";
  for (unsigned int i = 0; i < scores.size(); i++) {
    ScorePrint& score = scores[i];
    score_out << score.contig1 << "\t" << score.contig2 << "\t"
	      << score.score << endl;
  }
  score_out.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  clock_t start_t = clock();

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

  // init leaf nodes
  init_nodes(g_contigs, g_nodes);

  // set contig and node connectivity
  vector<ScorePrint> scores;
  process_contacts(contig_map, g_contigs, g_nodes, contacts, scores, params.stype, g_matrix);
  output_scores(params, scores);

  NodePairQueue queue;

  // fill queue with all connected nodes
  cout << "init queue..." << endl;
  init_queue(queue, params.min_score, params.correct_cluster_size);

  cout << "queue start size: " << queue.size() << endl;

  clock_t prev_t = clock();

  int count = g_nodes.size();
  while (queue.size() > 0) {
    // get first in queue
    NodePair top = queue.top();
    queue.pop();

    //    cout << "top i1=" << top.node_index1 << ", i2=" << top.node_index2 << ", s=" << top.score << endl;

    Node& node1 = get_node(top.node_index1);
    Node& node2 = get_node(top.node_index2);

    if (!node1.root_flag || !node2.root_flag)
      continue;

    double score = top.score;

    Node& new_node = add_node(top.score);
    merge_nodes(new_node, node1, node2, params.min_score, params.max_contig_pairs, params.stype, params.correct_cluster_size);

    // add new candidate pairs to queue
    int new_items = update_queue(queue, new_node, params.min_score, params.correct_cluster_size);

    clock_t current_t = clock();
    if ((((float)current_t-(float)prev_t) / CLOCKS_PER_SEC)> 600) {
      prev_t = current_t;
      double delta = (((float)current_t-(float)start_t) / CLOCKS_PER_SEC) / 60;
      cout << "N=" << count << ", min=" << delta << ", S=" << score << " ,Q=" << queue.size() << "(+" << new_items << ")" << endl;
    }

    count--;
  }

  // when done go over nodes and collect all roots
  NodeQueue root_queue;
  for (unsigned int i=0; i<g_nodes.size(); i++) {
    Node* node = &get_node(i);
    if (node->root_flag)
      root_queue.push(node);
  }

  // sorted roots
  vector<int> roots;
  while (root_queue.size() > 0) {
    roots.push_back(root_queue.top()->index);
    root_queue.pop();
  }

  cout << "number of trees when stopping: " << roots.size() << endl;

  output_results(params, roots);

  clock_t end_t = clock();
  double total_min = (((float)end_t-(float)start_t) / CLOCKS_PER_SEC) / 60;
  cout << "total run time (minutes): " << total_min << endl;

  return (0);
}
