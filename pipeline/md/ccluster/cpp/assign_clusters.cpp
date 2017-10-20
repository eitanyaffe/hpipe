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

struct Contig
{
  // string iden
  string name;

  // index
  int index;
};

// node in cluster tree
struct Node
{
  // index
  int index;

  // all contained contigs
  vector<string> contigs;

  // node children
  bool leaf ,root;
  int left, right;

  // nodes were united based on this score
  double score;
};

// dynamic list of cluster nodes
map<int, Node> g_nodes;

Node& get_node(unsigned int index)
{
  massert(g_nodes.find(index) != g_nodes.end(), "node not found");
  return (g_nodes[index]);
}

vector<string> set_contigs(Node& node, map<int, Node>& nodes, map<int, string>& contig_map)
{
  vector<string>& contigs = node.contigs;
  if (node.leaf) {
    massert(contig_map.find(node.index) != contig_map.end(), "contig of node not found");
    massert(contigs.size() == 0, "node should be empty");
    contigs.push_back(contig_map[node.index]);
  } else {
    Node& left = nodes[node.left];
    Node& right = nodes[node.right];

    vector<string> lc = set_contigs(left, nodes, contig_map);
    vector<string> rc = set_contigs(right, nodes, contig_map);

    contigs.insert(contigs.end(), lc.begin(), lc.end());
    contigs.insert(contigs.end(), rc.begin(), rc.end());
  }

  return contigs;
}

void set_clusters(Node& node, map<int, Node>& nodes, int min_elements, double min_score,
		  map<string, int>& contig_clusters, int& cluster, int level)
{
  if (node.leaf)
    return;

  Node& left = nodes[node.left];
  Node& right = nodes[node.right];

  if (((int)node.contigs.size() < min_elements) || (node.score/node.contigs.size() > min_score/100)) {
    //    cout << "cluster=" << cluster << ", level=" << level << ", node=" << node.index << " ,N=" << node.contigs.size() << ", score=" << node.score << ", lscore=" << left.score << ", rscore=" << right.score << endl;

    for (unsigned int i=0; i < node.contigs.size(); i++)
      contig_clusters[node.contigs[i]] = cluster;
    cluster++;
    return;
  }

  if ((int)left.contigs.size() > min_elements)
    set_clusters(left, nodes, min_elements, min_score, contig_clusters, cluster, level+1);

  if ((int)right.contigs.size() > min_elements)
    set_clusters(right, nodes, min_elements, min_score, contig_clusters, cluster, level+1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Parsing and handling user arguments
//////////////////////////////////////////////////////////////////////////////////////////////////

class UserParams {
public:
  string contig_ifn; // contig table
  string tree_ifn; // contact table

  string contig_ofn;

  double min_score;
  double min_elements;

  // ctor with default values
  UserParams() {
    contig_ifn = "c";
    tree_ifn = "t";
    min_score = 4;
    min_elements = 10;
    contig_ofn = "contigs.out";
  }
};

void usage(const char* name, UserParams& params)
{
  fprintf(stderr, "usage: %s [options]\n", name);
  cout << " -contigs <fn>: contig table" << endl;
  cout << " -tree <fn>: tree table" << endl;
  cout << " -min_score <double>: average percentage of shared neighbours over which we define a cell" << endl;
  cout << " -min_elements <int>: minimal number of contigs in cluster" << endl;
  cout << " -out_contigs <fn>: contig table output" << endl;
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
      else if (option == "-tree")
	params.tree_ifn = arg;
      else if (option == "-min_score")
	params.min_score = atof(arg);
      else if (option == "-out_contigs")
	params.contig_ofn = arg;
      else if (option == "-min_elements")
	params.min_elements = atoi(arg);
      else {
	cout << "Error: unknown option: " << option << endl;
	exit(1);
      }

      i += 2;
    }
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

void read_contigs(string fn, map<int, string>& contig_map, vector<string>& contigs)
{
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int id_idx = get_field_index("contig", fields);
  int node_idx = get_field_index("node", fields);
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string name = fields[id_idx];
    int node_index = atoi(fields[node_idx].c_str());

    // save map from name to index for contact table
    contig_map[node_index] = name;

    contigs.push_back(name);
  }
  in.close();
}

void read_tree(string fn, map<int, Node>& nodes, vector<int>& roots)
{
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int node_idx = get_field_index("node", fields);
  int level_idx = get_field_index("level", fields);
  int leaf_idx = get_field_index("leaf", fields);
  int score_idx = get_field_index("score", fields);
  int left_idx = get_field_index("left", fields);
  int right_idx = get_field_index("right", fields);
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    Node node;

    node.index = atoi(fields[node_idx].c_str());
    node.root = fields[level_idx] == "0";
    node.leaf = fields[leaf_idx] == "1";
    node.score = atof(fields[score_idx].c_str());
    node.left = atoi(fields[left_idx].c_str());
    node.right = atoi(fields[right_idx].c_str());

    nodes[node.index] = node;

    if (node.root)
      roots.push_back(node.index);
  }
  in.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  UserParams params;
  parse_user_arguments(argc, argv, params);

  // map from leaf node index to contig string id
  map<int, string> contig_map;
  vector<string> contigs;

  // init contigs and contacts
  cout << "contig table: " << params.contig_ifn << endl;
  read_contigs(params.contig_ifn, contig_map, contigs);

  // roots are indices into g_nodes
  vector<int> roots;

  cout << "tree table: " << params.tree_ifn << endl;
  read_tree(params.tree_ifn, g_nodes, roots);

  for (unsigned int i=0; i<roots.size(); i++) {
    Node& node = get_node(roots[i]);
    set_contigs(node, g_nodes, contig_map);
  }

  // map from contig to cluster
  map<string, int> contig_clusters;
  int cluster = 1;

  for (unsigned int i=0; i<roots.size(); i++) {
    Node& node = get_node(roots[i]);
    if ((int)node.contigs.size() > params.min_elements)
      set_clusters(node, g_nodes, params.min_elements, params.min_score, contig_clusters, cluster, 0);
  }
  cout << "number of clusters: " << cluster-1 << endl;

  // finally output
  cout << "writing result: " << params.contig_ofn << endl;
  ofstream out(params.contig_ofn.c_str());
  massert(out.is_open(), "could not open file %s", params.contig_ofn.c_str());
  out << "contig\tcluster\n";
  for (unsigned int i = 0; i < contigs.size(); i++) {
    string name = contigs[i];
    int cluster = (contig_clusters.find(name) != contig_clusters.end()) ?  contig_clusters[name] : -1;
    out << name << "\t" << cluster << endl;
  }
  out.close();

  return (0);
}
