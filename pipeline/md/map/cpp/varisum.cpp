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

enum VariType { vtNone, vtSubstitute, vtDelete, vtInsert, vtDangleLeft, vtDangleRight };
struct Variation {
  VariType type;

  // vtSubstitute (N=1), vtInsert (N>=1)
  string seq;

  // ctor
  Variation(VariType _type=vtNone, string _seq="NA") : type(_type), seq(_seq) {};
  bool get_type() { return type; };
  string type_str() {
    switch(type) {
    case vtNone: return "none";
    case vtSubstitute: return "snp";
    case vtDelete: return "delete";
    case vtInsert: return "insert";
    case vtDangleLeft: return "dangle_left";
    case vtDangleRight: return "dangle_right";
    }
    return "";
  }
  string str() {
    return type_str() + "_" + seq;
  }
};

bool operator<(const Variation& lhs, const Variation& rhs) {
  if (lhs.type != rhs.type)
    return (int)lhs.type < (int)rhs.type;
  return lhs.seq < rhs.seq;

}

void init_params(int argc, char **argv, Parameters& params)
{
  params.add_parser("idir", new ParserFilename("input directory"), true);
  params.add_parser("contigs", new ParserFilename("contig table"), true);
  params.add_parser("margin", new ParserInteger("safety margin from read edges (nts)"), true);
  params.add_parser("ofn_nt", new ParserFilename("output nt table"), true);
  params.add_parser("ofn_link", new ParserFilename("output nt link table"), true);

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

void parse_cigar(string in, vector<pair < char, int> > & result)
{
  result.resize(0);
  string field;
  for (unsigned int i=0; i<in.size(); ++i) {
    char c = in[i];
    int length = 0;
    if (isdigit(c))
      length = length*10 + (c - '0');
    else
      result.push_back(make_pair(c, length));
  }
}

string reverse_complement(string seq) {
  string result = seq;
  int N = seq.length();
  for (int i=0; i<N; i++) {
    char c = seq[N-i-1], r;
    switch( c )
      {
      case 'A': r = 'T'; break;
      case 'G': r = 'C'; break;
      case 'C': r = 'G'; break;
      case 'T': r = 'A'; break;
      default: r = c;
      }
    result[i] = r;
  }
  return(result);
}

void process_fasta(string fn,
		   map<string, int>& contig_map,
		   int margin,
		   map< string, map< int, map <Variation, int> > >& nt_coord,
		   map< string, map< int, map <Variation, int> > >& link_coord,
		   map<string, vector<int> >& total_nt,
		   map<string, vector<int> >& total_link)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  vector<string> sub_fields;
  char sub_delim = ';';

  vector<string> ssub_fields;
  char ssub_delim = ',';

  // parse header
  split_line(in, fields, delim);
  int id_ind = get_field_index("id", fields);
  int score_ind = get_field_index("score", fields);
  int contig_ind = get_field_index("contig", fields);
  int coord_ind = get_field_index("coord", fields);
  int back_coord_ind = get_field_index("back_coord", fields);
  int strand_ind = get_field_index("strand", fields);
  int sub_ind = get_field_index("substitute", fields);
  int insert_ind = get_field_index("insert", fields);
  int delete_ind = get_field_index("delete", fields);
  int cigar_ind = get_field_index("cigar", fields);

  map<string, int> id_score;

  int read_count = 1;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    if (read_count++ % 1000000 == 0) {
      cout << "read: " << read_count-1 << " ..." << endl;
    }

    // use only longest match for each read
    string id = fields[id_ind];

    int score = atoi(fields[score_ind].c_str());
    if (id_score.find(id) != id_score.end()) {
      massert(score <= id_score[id], "first line per id should be the best score");
      continue;
    }
    id_score[id] = score;

    string contig = fields[contig_ind];
    int front_coord = atoi(fields[coord_ind].c_str()) - 1;
    int back_coord = atoi(fields[back_coord_ind].c_str()) - 1;
    bool strand = fields[strand_ind] == "1";
    int left_coord = strand ? back_coord : front_coord;
    int right_coord = !strand ? back_coord : front_coord;
    int read_length = right_coord - left_coord + 1;

    // skip short contigs which are not in contig table
    if(contig_map.find(contig) == contig_map.end())
      continue;

    int contig_length = contig_map[contig];

    string sub = fields[sub_ind];
    string insert = fields[insert_ind];
    string del = fields[delete_ind];

    // handle left/right dangle
    string cigar_str = fields[cigar_ind];
    vector<pair < char, int> > cigar;
    parse_cigar(cigar_str, cigar);

    map< int, map <Variation, int> >& nt_coord_contig = nt_coord[contig];
    map< int, map <Variation, int> >& link_coord_contig = link_coord[contig];

    bool clip_start = (cigar.front().first == 'S' || cigar.front().first == 'H') && (left_coord != 0);
    bool clip_end = (cigar.back().first == 'S' || cigar.back().first == 'H') && (right_coord != contig_length-1);
    //bool clip_start = (cigar.front().first == 'S' || cigar.front().first == 'H');
    //bool clip_end = (cigar.back().first == 'S' || cigar.back().first == 'H');

    vector<int>& total_nt_contig = total_nt[contig];
    vector<int>& total_link_contig = total_link[contig];

    //cout << "id=" << id << ", front_coord=" << front_coord << ", back_coord=" << back_coord << ", strand=" << (strand ? "+" : "-") << endl;
    // cout << "id=" << id << ", left_coord=" << left_coord << ", right_coord=" << right_coord << ", length=" << read_length << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del << endl;

    // dangle left
    if (strand ? clip_start : clip_end) {
      int coord = strand ? left_coord : right_coord;
      Variation var(strand ? vtDangleLeft : vtDangleRight);
      INCREMENT(link_coord_contig, coord, var, contig_length);
      // cout << "id=" << id << ", contig=" << contig << ", left_coord=" << left_coord << ", right_coord=" << right_coord
      // 	   << ", strand=" << (strand ? "+" : "-") << ", length=" << read_length
      // 	   << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del << endl;
    }

    // dangle right
    if (!strand ? clip_start : clip_end) {
      int coord = strand ? right_coord : left_coord;
      Variation var(strand ? vtDangleRight : vtDangleLeft);
      INCREMENT(link_coord_contig, coord, var, contig_length);
      // cout << "id=" << id << ", contig=" << contig << ", left_coord=" << left_coord << ", right_coord=" << right_coord
      // 	   << ", strand=" << (strand ? "+" : "-") << ", length=" << read_length
      // 	   << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del << endl;
    }

    // substitutes
    split_string(sub, sub_fields, sub_delim);
    for (unsigned int i=0; i<sub_fields.size(); ++i) {
      split_string(sub_fields[i], ssub_fields, ssub_delim);
      if (ssub_fields.size() < 4)
	continue;
      int coord = atoi(ssub_fields[0].c_str()) - 1;
      if (coord < left_coord+margin || coord > right_coord-margin)
	continue;
      // string nt = (strand ? ssub_fields[3] : reverse_complement(ssub_fields[3]));
      string nt = ssub_fields[3];
      if (nt == "N")
	continue;
      Variation var(vtSubstitute, nt);
      INCREMENT(nt_coord_contig, coord, var, contig_length);
      // cout << "id=" << id << ", contig=" << contig << ", left_coord=" << left_coord << ", right_coord=" << right_coord
      // 	   << ", strand=" << (strand ? "+" : "-") << ", length=" << read_length
      // 	   << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del << endl;
    }

    // deletions
    split_string(del, sub_fields, sub_delim);
    for (unsigned int i=0; i<sub_fields.size(); ++i) {
      split_string(sub_fields[i], ssub_fields, ssub_delim);
      if (ssub_fields.size() < 2)
	continue;
      int start_coord = atoi(ssub_fields[0].c_str()) - 1;
      if (start_coord < left_coord+margin || start_coord > right_coord-margin)
	continue;
      int del_length = atoi(ssub_fields[1].c_str());
      for (int k=0; k<del_length; k++) {
	int coord = start_coord + k;
	Variation var(vtDelete);
	INCREMENT(nt_coord_contig, coord, var, contig_length);
      }
      // cout << "id=" << id << ", contig=" << contig << ", left_coord=" << left_coord << ", right_coord=" << right_coord
      // 	   << ", strand=" << (strand ? "+" : "-") << ", length=" << read_length
      // 	   << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del << endl;
    }

    // insertions
    split_string(insert, sub_fields, sub_delim);
    for (unsigned int i=0; i<sub_fields.size(); ++i) {
      split_string(sub_fields[i], ssub_fields, ssub_delim);
      if (ssub_fields.size() < 3)
	continue;
      int coord = atoi(ssub_fields[0].c_str()) - 2;
      if (coord < left_coord+margin || coord > right_coord-margin)
	continue;
      // string seq = (strand ? ssub_fields[2] : reverse_complement(ssub_fields[2]));
      string seq = ssub_fields[2];
      Variation var(vtInsert, seq);
      INCREMENT(link_coord_contig, coord, var, contig_length);
      // cout << "id=" << id << ", contig=" << contig << ", left_coord=" << left_coord << ", right_coord=" << right_coord
      // 	   << ", strand=" << (strand ? "+" : "-") << ", length=" << read_length
      // 	   << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del << endl;
    }

    for (int i=margin; i<read_length-margin; ++i) {
      int coord = left_coord + i;
      massert_range(total_nt_contig, coord);
      total_nt_contig[coord]++;
      if (i+1 < read_length)
	total_link_contig[coord]++;
    }
  }
}

void read_contig_table(string fn, map<string, int>& contig_map)
{
  cout << "reading table: " << fn << endl;
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
    contig_map[contig] = length;
  }
}

inline double mround(double v, int digits) {
  double f = pow(10,digits);
  return (round(v * f) / f);
}

void save_table(map<string, int>& contig_map,
		map< string, map< int, map <Variation, int> > >& table,
		map<string, vector<int> >& total,
		string ofn)
{
  cout << "writing output to table: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  out << "contig\tcoord\ttype\tcount\tpercent\tsequence" << endl;

  for (map<string, int>::iterator it=contig_map.begin(); it!=contig_map.end(); ++it) {
    string contig = (*it).first;
    map< int, map <Variation, int> >& table_contig = table[contig];
    vector<int>& total_contig = total[contig];

    for (map< int, map <Variation, int> >::iterator it = table_contig.begin(); it != table_contig.end(); ++it) {
      int coord = (*it).first;
      map <Variation, int>& xmap = (*it).second;
      int vcount = 0;
      for (map <Variation, int>::iterator jt=xmap.begin(); jt!=xmap.end(); ++jt) {
	Variation var = (*jt).first;
	int count =  (*jt).second;
	vcount += count;
	out << contig << "\t" << coord+1 << "\t" << var.type_str() << "\t" << count << "\t" << mround(100*count/(double)total_contig[coord],1) << "\t" << var.seq << endl;
      }
      int rcount = (total_contig[coord] - vcount);
      out << contig << "\t" << coord+1 << "\tREF\t" << rcount << "\t" << mround(100*rcount/(double)total_contig[coord],1) << "\tNA" << endl;
    }
  }
  out.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  Parameters params;
  init_params(argc, argv, params);

  string contig_table_ifn = params.get_string("contigs");
  int margin = params.get_int("margin");
  string idir = params.get_string("idir");
  string ofn_nt = params.get_string("ofn_nt");
  string ofn_link = params.get_string("ofn_link");

  map<string, int> contig_map;
  read_contig_table(contig_table_ifn, contig_map);
  cout << "number of contigs: " << contig_map.size() << endl;

  // for deletions and substitutes
  map< string, map< int, map <Variation, int> > > nt_coord;
  map<string, vector<int> > total_nt;

  // for inserts and dangles
  map< string, map< int, map <Variation, int> > > link_coord;
  map<string, vector<int> > total_link;

  for (map<string, int>::iterator it=contig_map.begin(); it!=contig_map.end(); ++it) {
    string contig = (*it).first;
    int length = (*it).second;

    vector<int>& total_nt_contig = total_nt[contig];
    vector<int>& total_link_contig = total_link[contig];
    total_nt_contig.resize(length);
    total_link_contig.resize(length);
  }

  DIR *dir;
  struct dirent *ent;
  massert((dir = opendir (idir.c_str())) != NULL, "could not open input directory");
  // int file_count = 0;
  while ((ent = readdir (dir)) != NULL) {
    string ifn = ent->d_name;
    if (ifn.find("fast") == string::npos || ifn.find("~") != string::npos)
      continue;
    process_fasta(idir + "/" + ifn, contig_map, margin, nt_coord, link_coord, total_nt, total_link);
    //    if (file_count++ > 3)
    //      break;
  }
  closedir (dir);

  save_table(contig_map, nt_coord, total_nt, ofn_nt);
  save_table(contig_map, link_coord, total_link, ofn_link);
}
