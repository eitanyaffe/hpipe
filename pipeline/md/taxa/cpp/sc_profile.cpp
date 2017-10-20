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

#include <dirent.h>

#include "util.h"
#include "Params.h"

#define massert_range(v, coord) massert(coord >= 0 && coord < (int)v.size(), "coord not in vector range");
// #define INCREMENT(v, coord, var, length) massert(coord >=0 && coord < length, "coord out of range, coord=%d, contig_size=%d", coord, length); v[coord][var]++; cout << coord+1 <<  " : " << var.str() << endl;
#define INCREMENT(v, coord, var, length) massert(coord >=0 && coord < length, "coord out of range, coord=%d, contig_size=%d", coord, length); v[coord][var]++;

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
// Coverage class
///////////////////////////////////////////////////////////////////////////////////////

class Coverage {
private:
  int m_length;
  map< string, int > m_contig_length;
  map< string, vector<int> > m_contig_edit;
public:
  Coverage(int length) : m_length(length) {};
  void read_contig_table(string fn);
  void init(string fn);
  bool is_valid(string contig, int start, int end);
  void append(string contig, int start, int end, int edit);
  void save_table(string ofn);
  void save_summary(string ofn);
};

void Coverage::read_contig_table(string fn)
{
  cout << "reading contig table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int contig_ind = get_field_index("contig", fields);
  int length_ind = get_field_index("length", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string contig = fields[contig_ind];
    int length = atoi(fields[length_ind].c_str());
    m_contig_length[contig] = length;
  }
}

void Coverage::init(string fn)
{
  read_contig_table(fn);
  for (map<string, int>::iterator it=m_contig_length.begin(); it!=m_contig_length.end(); ++it) {
    string contig = (*it).first;
    int length = (*it).second;
    vector<int>& edit = m_contig_edit[contig];
    edit.resize(length);
    for (int i=0; i<length; ++i)
      edit[i] = -1;
  }
}

void Coverage::save_table(string ofn)
{
  cout << "writing summary table: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  out << "contig\tcoord\tedit_distance" << endl;

  for (map<string, int>::iterator it=m_contig_length.begin(); it!=m_contig_length.end(); ++it) {
    string contig = (*it).first;
    int contig_length = (*it).second;
    vector < int > & edit_v = m_contig_edit[contig];
    for (int coord = 0; coord < contig_length; coord++)
      out << contig << "\t" << coord+1 << "\t" << edit_v[coord] << endl;
  }
  out.close();
}

void Coverage::save_summary(string ofn)
{
  map<int, int> result;
  for (map<string, int>::iterator it=m_contig_length.begin(); it!=m_contig_length.end(); ++it) {
    string contig = (*it).first;
    int contig_length = (*it).second;
    vector < int > & edit_v = m_contig_edit[contig];
    for (int coord = 0; coord < contig_length; coord++)
      result[edit_v[coord]]++;
  }

  cout << "writing output table: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  out << "edit\tcount" << endl;
  for (map<int, int>::iterator it=result.begin(); it!=result.end(); ++it) {
    int edit = (*it).first;
    int count = (*it).second;
    out << edit << "\t" << count << endl;
  }

  out.close();
}

bool Coverage::is_valid(string contig, int start, int end)
{
  int length = end - start + 1;
  return m_contig_length.find(contig) != m_contig_length.end() && length == m_length;
}

void Coverage::append(string contig, int start, int end, int edit)
{
  int length = end - start + 1;
  if (m_contig_length.find(contig) == m_contig_length.end() || length != m_length)
    return;

  int contig_length = m_contig_length[contig];
  massert(start >=0 && end < contig_length, "appended read out of range, contig=%s start=%d, end=%d, length=%d", contig.c_str(), start, end, contig_length);
  vector<int>& edit_v = m_contig_edit[contig];

  // attribute segment to start coord
  for (int coord=start; coord<=end; ++coord)
    if (edit < edit_v[coord] || edit_v[coord] == -1)
      edit_v[coord] = edit;
}

///////////////////////////////////////////////////////////////////////////////////////

void init_params(int argc, char **argv, Parameters& params)
{
  params.add_parser("idir", new ParserFilename("input directory"), true);
  params.add_parser("parse_src", new ParserBoolean("read id includes source contig/coord"), true);
  params.add_parser("src_contig_table", new ParserFilename("source contig table"), false);
  params.add_parser("tgt_contig_table", new ParserFilename("target contig table"), true);
  params.add_parser("read_length", new ParserInteger("read length"), true);
  params.add_parser("src_ofn", new ParserFilename("src output table"), false);
  params.add_parser("tgt_ofn", new ParserFilename("tgt output table"), true);
  params.add_parser("src_summary_ofn", new ParserFilename("src output summary table"), false);
  params.add_parser("tgt_summary_ofn", new ParserFilename("tgt output summary table"), true);

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

void split_string(string in, vector<string> &fields, char delim)
{
  fields.resize(0);
  string field;
  for (unsigned int i=0; i<in.size(); ++i) {
    char c = in[i];
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

inline int parse_coord(string str)
{
  return atoi(str.substr(str.find(':')+1, string::npos).c_str());
}

void process_sam(string fn, bool parse_src, Coverage* src_coverage, Coverage* tgt_coverage)
{
  cout << "processing sam table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);

  int src_id_ind = get_field_index("id", fields);
  int tgt_contig_ind = get_field_index("contig", fields);
  int tgt_coord_ind = get_field_index("coord", fields);
  int tgt_back_coord_ind = get_field_index("back_coord", fields);
  int tgt_strand_ind = get_field_index("strand", fields);
  int edit_ind = get_field_index("edit_dist", fields);

  int read_count = 1;

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    if (read_count++ % 1000000 == 0)
      cout << "read: " << read_count-1 << " ..." << endl;

    string src_id = fields[src_id_ind];

    string tgt_contig = fields[tgt_contig_ind];
    int tgt_front_coord = atoi(fields[tgt_coord_ind].c_str()) - 1;
    int tgt_back_coord = atoi(fields[tgt_back_coord_ind].c_str()) - 1;
    bool tgt_strand = fields[tgt_strand_ind] == "1";
    int edit = atoi(fields[edit_ind].c_str());

    int tgt_start = tgt_strand ? tgt_back_coord : tgt_front_coord;
    int tgt_end = !tgt_strand ? tgt_back_coord : tgt_front_coord;

    if(!tgt_coverage->is_valid(tgt_contig, tgt_start, tgt_end))
      continue;

    if (parse_src) {
      vector<string> id_fields;
      split_string(fields[src_id_ind], id_fields, '_');
      massert(id_fields.size() == 3, "source id must have 3 fields: %s", fields[src_id_ind].c_str());
      string src_contig = id_fields[0];
      int src_start = parse_coord(id_fields[1]) - 1;
      int src_end = parse_coord(id_fields[2]) - 1;
      if(!src_coverage->is_valid(src_contig, src_start, src_end))
	continue;
      src_coverage->append(src_contig, src_start, src_end, edit);
      tgt_coverage->append(tgt_contig, tgt_start, tgt_end, edit);
    } else {
      tgt_coverage->append(tgt_contig, tgt_start, tgt_end, edit);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  Parameters params;
  init_params(argc, argv, params);

  string idir = params.get_string("idir");
  bool parse_src = params.get_bool("parse_src");

  string tgt_contig_table_ifn = params.get_string("tgt_contig_table");
  int read_length = params.get_int("read_length");
  string tgt_ofn = params.get_string("tgt_ofn");
  string tgt_summary_ofn = params.get_string("tgt_summary_ofn");

  Coverage *tgt_coverage, *src_coverage;

  tgt_coverage = new Coverage(read_length);
  tgt_coverage->init(tgt_contig_table_ifn);

  cout << "searching for sam files in directory: " << idir << endl;

  if (parse_src) {
    src_coverage = new Coverage(read_length);
    string src_contig_table_ifn = params.get_string("src_contig_table");
    src_coverage->init(src_contig_table_ifn);
  } else {
    src_coverage = NULL;
  }

  DIR *dir;
  struct dirent *ent;
  massert((dir = opendir (idir.c_str())) != NULL, "could not open input directory");

  while ((ent = readdir (dir)) != NULL) {
    string ifn = ent->d_name;
    if (ifn == "." || ifn == "..")
      continue;
    process_sam(idir + "/" + ifn, parse_src, src_coverage, tgt_coverage);
  }
  closedir (dir);

  tgt_coverage->save_table(tgt_ofn);
  tgt_coverage->save_summary(tgt_summary_ofn);

  if (parse_src) {
    string src_summary_ofn = params.get_string("src_summary_ofn");
    string src_ofn = params.get_string("src_ofn");
    src_coverage->save_table(src_ofn);
    src_coverage->save_summary(src_summary_ofn);
  }
}
