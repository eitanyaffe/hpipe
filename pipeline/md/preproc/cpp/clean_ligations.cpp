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

#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>

#include <ctime>
#include <unordered_map>

#include "util.h"
#include "Kmer.h"

using namespace std;

struct FragmentEnd
{
  // ligation neighbour counts
  map<Kmer, int> neighbor_map;

  // post process fields
  bool has_genomic_neighbor;
  Kmer genomic_neighbor_kmer;

  FragmentEnd() :  has_genomic_neighbor(false) {};
};

typedef unordered_map<Kmer, FragmentEnd, KHash<size_t> > LigationMap;

enum Mode { mSplit, mTrim };

string mode2str(Mode mode)
{
  switch(mode) {
  case mSplit: return("split");
  case mTrim: return("trim");
  default: return("unknown");
  }
}

struct Stats
{
  // reads
  int read_count;

  // sites
  int site_count;
  int site_genomic_count;

  // segments
  int segment_count;
  int short_segment_count;

  // fends
  int fend_count;
  int fend_genomic_count;

  Stats() : read_count(0), site_count(0), site_genomic_count(0), segment_count(0), short_segment_count(0), fend_count(0), fend_genomic_count(0) {};
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Parsing and handling user arguments
//////////////////////////////////////////////////////////////////////////////////////////////////

class UserParams {
public:
  Mode mode;
  int ksize;
  string site;
  int min_coverage;
  int min_read_length;
  double min_ratio;
  string ofn_stats;
  
  string input_dir;
  string side1, side2;
  string output_dir;
  UserParams() : mode(mSplit), ksize(96), site("AGCT"), min_coverage(3), min_read_length(50), min_ratio(2), ofn_stats("ostats"), side1("R1"), side2("R2")  {};
};

void usage(const char* name, UserParams& params)
{
  fprintf(stderr, "clean_ligation: Identify 3C ligation junctions. Handle junctions according to mode\n");
  fprintf(stderr, "usage: %s [options]\n", name);
  cout << " -mode <strict|split|single>: work mode" << endl;
  cout << "   split: breakdown read by ligation junctions, disrupting read pairing" << endl;
  cout << "   trim: trim up to first ligation junction on read, maintaining valid read pairing" << endl;
  cout << " -input_dir<path>: read input dir" << endl;
  cout << " -input_side1<string>: substring indicator in filename of side 1" << endl;
  cout << " -input_side2<string>: substring indicator in filename of side 2" << endl;
  cout << " -output_dir<path>: read output dir" << endl;
  cout << " -ksize <int>: k size" << endl;
  cout << " -site <string>: restriction cut site" << endl;
  cout << " -min_coverage <int>: minimal coverage of ligation junction" << endl;
  cout << " -min_ratio <double>: minimal enrichment over best competing junction" << endl;
  cout << " -min_read_length <int>: threshold on output read length, relevant only for split mode" << endl;
  cout << " -ofn_stats <fn>: output stats" << endl;
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

      if (option == "-input_dir")
	params.input_dir = arg;
      else if (option == "-input_side1")
	params.side1 = arg;
      else if (option == "-input_side2")
	params.side2 = arg;
      else if (option == "-output_dir")
	params.output_dir = arg;
      else if (option == "-ksize")
	params.ksize = atoi(arg);
      else if (option == "-mode") {
	string mode = string(arg);
	if (mode == "trim")
	  params.mode = mTrim;
	else if (mode == "split")
	  params.mode = mSplit;
	else {
	  cout << "Error: unknown mode: " << mode << endl;
	  exit(1);
	}
      } else if (option == "-min_read_length")
	params.min_read_length = atoi(arg);
      else if (option == "-site")
	params.site = string(arg);
      else if (option == "-min_coverage")
	params.min_coverage = atoi(arg);
      else if (option == "-min_ratio")
	params.min_ratio = atof(arg);
      else if (option == "-ofn_stats")
	params.ofn_stats = arg;
      else {
	cout << "Error: unknown option: " << option << endl;
	exit(1);
      }

      i += 2;
    }
}

ifstream::pos_type filesize(string filename)
{
    ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

void populate_ligation_map(ifstream& in, LigationMap& lmap, string& site, unsigned int ksize)
{
  const int MAXLINE = 1024;
  long counter = 0;
  char line[4*MAXLINE];

  int site_length = site.length();
  while(1) {
    int index = counter % 4;
    in.getline(line + index*MAXLINE, MAXLINE-1);
    bool eof = in.eof();

    if (!eof && !in.good()) {
      cerr << "Error reading line: " << counter+1 << ", rdstate=" << in.rdstate() << endl;
      if ( (in.rdstate() & std::ifstream::failbit ) != 0)
	cerr << "Probably line in file exceeds allocated maximal length: " << MAXLINE << endl;
      exit(-1);
    }

    // check if sequence already encountered
    if (index == 1) {
      string seq(line + MAXLINE);
      // cout << seq << endl;
      for (unsigned int i=0; i<seq.length()-site_length; i++) {
	// cout << "i=" << i << ", len=" << seq.length() << ", sub_len=" << site_length << ", ksize=" << ksize << endl;
	string psite = seq.substr(i, site_length);
	if (!(psite == site) || (i < ksize) || ((i+site_length+ksize) > seq.length()))
	  continue;

	string fwd = seq.substr(i-ksize, ksize);
	string bck = seq.substr(i+site_length, ksize);
	if (fwd.find('N') != string::npos || bck.find('N') != string::npos)
	  continue;
	Kmer fwd_kmer = Kmer(fwd);
	Kmer bck_kmer = Kmer(reverse_complement(bck));

	lmap[fwd_kmer].neighbor_map[bck_kmer]++;
	lmap[bck_kmer].neighbor_map[fwd_kmer]++;
      }
    }

    if (eof)
      break;

    counter++;
  }
}

void detach_ligations(ifstream& in, ofstream& out, LigationMap& lmap, string& site, unsigned int ksize, Mode mode, int min_read_length, string side, Stats& stats)
{
  const int MAXLINE = 1024;
  long counter = 0;
  char line[4*MAXLINE];
  int site_length = site.length();

  while(1) {
    int index = counter % 4;
    in.getline(line + index*MAXLINE, MAXLINE-1);
    bool eof = in.eof();

    if (!eof && !in.good()) {
      cerr << "Error reading line: " << counter+1 << ", rdstate=" << in.rdstate() << endl;
      if ( (in.rdstate() & std::ifstream::failbit ) != 0)
	cerr << "Probably line in file exceeds allocated maximal length: " << MAXLINE << endl;
      exit(-1);
    }

    // check if sequence already encountered
    if (index == 3) {
      string seq(line + MAXLINE);
      string quality(line + 3*MAXLINE);
      vector<int> detach_coords;

      stats.read_count++;

      for (unsigned int i=0; i<seq.length(); i++) {
	if (!(seq.substr(i, site_length) == site))
	  continue;
	stats.site_count++;

	if ((ksize <= i) && ((i+site_length+ksize) <= seq.length())) {
	  
	  string fwd = seq.substr(i-ksize, ksize);
	  string bck = seq.substr(i+site_length, ksize);
	  if (fwd.find('N') != string::npos || bck.find('N') != string::npos)
	    continue;
	  Kmer fwd_kmer = Kmer(seq.substr(i-ksize, ksize));
	  Kmer bck_kmer = Kmer(reverse_complement(seq.substr(i+site_length, ksize)));

	  // skip if site is considered bona fide genomic
	  if (!(lmap.find(fwd_kmer) != lmap.end() &&
		lmap.find(bck_kmer) != lmap.end() &&
		lmap[fwd_kmer].has_genomic_neighbor &&
		lmap[bck_kmer].has_genomic_neighbor &&
		lmap[fwd_kmer].genomic_neighbor_kmer == bck_kmer &&
		lmap[bck_kmer].genomic_neighbor_kmer == fwd_kmer)) {
	    stats.site_genomic_count++;
	    continue;
	  }
	}
	detach_coords.push_back(i);
      }
      // output read, breakdowning down if needed
      if (detach_coords.size() == 0) {
	if (mode == mTrim)
	  out << line << endl;
	else
	  out << line << "_" << side << endl;
	for (int i=1; i<4; i++)
	  out << line + i*MAXLINE << endl;
      }	else {

	stats.segment_count += detach_coords.size() + 1;

	int prev_coord = 0;
	detach_coords.push_back(seq.length());
	for (unsigned int i=0; i<detach_coords.size(); i++) {
	  int coord = detach_coords[i];
	  int slength = coord - prev_coord;

	  if (slength < min_read_length || slength == 0) {
	    stats.short_segment_count++;
	    continue;
	  }

	  if (mode == mTrim)
	    out << string(line) << endl;
	  else
	    out << string(line) << "_" << side << "_" << i << endl;

	  out << seq.substr(prev_coord, slength) << endl;
	  out << line + 2*MAXLINE << endl;
	  out << quality.substr(prev_coord, slength) << endl;

	  prev_coord = coord + site_length;

	  if (mode == mTrim)
	    break;
	}
      }
    }
    if (eof)
      break;

    counter++;
  }
}

void get_paired_fns(string ipath, string opath,
		    string side1, string side2, 
		    vector<string>& ifns1, vector<string>& ifns2,
		    vector<string>& ofns1, vector<string>& ofns2)
{
  DIR *dir = opendir(ipath.c_str());
  massert(dir != NULL, "cannot locate input dir: %s", ipath.c_str());
  dirent *ent;

  set<string> keys;
  while ((ent = readdir (dir)) != NULL) {
    if (ent->d_type != DT_REG)
      continue;
    
    string fn = ent->d_name;
    if (fn.find(side1) != string::npos) {
      fn.replace(fn.find(side1), side1.length(), "^");
      keys.insert(fn);
    } else if (fn.find(side2) != string::npos) {
      fn.replace(fn.find(side2), side2.length(), "^");
      keys.insert(fn);
    }
  }
  closedir (dir);

  for (set<string>::iterator it = keys.begin(); it != keys.end(); ++it) {
    string fn1 = *it;
    string fn2 = *it;
    fn1.replace(fn1.find("^"), 1, side1);
    fn2.replace(fn2.find("^"), 1, side2);
    ifns1.push_back(ipath + "/" + fn1);
    ifns2.push_back(ipath + "/" + fn2);
    ofns1.push_back(opath + "/" + fn1);
    ofns2.push_back(opath + "/" + fn2);
  }

}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  UserParams params;
  LigationMap lmap;
  Stats stats;

  parse_user_arguments(argc, argv, params);

  massert(params.site != "", "site not defined");

  cout << "site: " << params.site << endl;
  cout << "ksize: " << params.ksize << endl;
  cout << "mode: " << mode2str(params.mode) << endl;
  cout << "min_coverage: " << params.min_coverage << endl;
  cout << "min_ratio: " << params.min_ratio << endl;
  cout << "min_read_length: " << params.min_read_length << endl;
  cout << "stats ofn: " << params.ofn_stats << endl;

  massert(params.ksize> 0, "ksize must be >0");
  massert(params.site == reverse_complement(params.site), "site not palindrome: %s", params.site.c_str());

  Kmer::set_k(params.ksize);

  vector<string> ifns1, ifns2, ofns1, ofns2;  
  get_paired_fns(params.input_dir, params.output_dir, 
		 params.side1, params.side2, 
		 ifns1, ifns2, ofns1, ofns2);
  massert(ifns1.size() == ifns2.size(), "side1 and side2 files do not match in counts");

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // populate ligation map
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "loading paired input files: " << ifns1.size() << endl;
  for (unsigned int i=0; i<ifns1.size(); i++) {
    string ifn1 = ifns1[i];
    string ifn2 = ifns2[i];

    cout << "reading file: " << ifn1 << endl;
    ifstream in1(ifn1.c_str());
    massert(in1.is_open(), "could not open file %s", ifn1.c_str());
    populate_ligation_map(in1, lmap, params.site, params.ksize);
    in1.close();

    cout << "reading file: " << ifn2 << endl;
    ifstream in2(ifn2.c_str());
    massert(in2.is_open(), "could not open file %s", ifn2.c_str());
    populate_ligation_map(in2, lmap, params.site, params.ksize);
    in2.close();

    cout << "."; cout.flush();
  }
  cout << endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // identify genomic junctions
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "all input files read into map, identifing ligation junctions..." << endl;
  for (LigationMap::iterator it = lmap.begin(); it != lmap.end(); ++it) {
    // Kmer kmer = (*it).first;
    FragmentEnd& fend = (*it).second;
    stats.fend_count++;

    // sort counts
    vector<int> counts;
    for (map<Kmer,int>::const_iterator jt = fend.neighbor_map.begin(); jt != fend.neighbor_map.end(); ++jt)
      counts.push_back((*jt).second);
    sort(counts.begin(), counts.end());
    int N = counts.size();
    int max_count = counts[N-1];
    int second_count = counts[N-2];
    // cout << "kmer: " << kmer.get_sequence() << ": n=" << fend.neighbor_map.size() << ", 1st=" << max_count << ", 2nd=" << second_count << endl;

    // skip if below thresholds
    if (max_count < params.min_coverage || (second_count > 0 && max_count/second_count < params.min_ratio)) {
      fend.has_genomic_neighbor = false;
      // cout << "non_genomic: " << kmer.get_sequence() << ", 1st: " << max_count << ", 2nd: " << second_count << ", ratio: " << max_count/second_count << endl;
      continue;
    }

    // find and mark putative genomic neighbor
    Kmer max_kmer;
    for (map<Kmer,int>::const_iterator jt = fend.neighbor_map.begin(); jt != fend.neighbor_map.end(); ++jt) {
      const Kmer& nkmer = (*jt).first;
      int count = (*jt).second;
      if (max_count == count) {
	max_kmer = nkmer;
	break;
      }
    }
    fend.has_genomic_neighbor = true;
    fend.genomic_neighbor_kmer = max_kmer;
    stats.fend_genomic_count++;

    // cout << "genomic neighbor found: " << max_kmer.get_sequence() << endl; 
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // clean ligations using map
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  cout << "second round: going over input files " << endl;
  for (unsigned int i=0; i<ifns1.size(); i++) {
    string ifn1 = ifns1[i];
    string ifn2 = ifns2[i];

    string ofn1 = ofns1[i];
    string ofn2 = ofns2[i];

    cout << "read1: " << ifn1 << " -> " << ofn1 << endl;
    cout << "read2: " << ifn2 << " -> " << ofn2 << endl;

    ofstream out1(ofn1.c_str());
    massert(out1.is_open(), "could not open file %s", ofn1.c_str());

    ofstream out2(ofn2.c_str());
    massert(out2.is_open(), "could not open file %s", ofn2.c_str());

    // cout << "input read1: " << ifn1 << endl;
    // cout << "input read2: " << ifn2 << endl;

    ifstream in1(ifn1.c_str());
    massert(in1.is_open(), "could not open file %s", ifn1.c_str());
    detach_ligations(in1, out1, lmap, params.site, params.ksize, params.mode, params.min_read_length, "R1", stats);
    in1.close();

    ifstream in2(ifn2.c_str());
    massert(in2.is_open(), "could not open file %s", ifn2.c_str());
    detach_ligations(in2, out2, lmap, params.site, params.ksize, params.mode, params.min_read_length, "R2", stats);
    in2.close();

    out1.close();
    out2.close();

    cout << "."; cout.flush();
  }
  cout << endl;
  
  cout << "output stats: " << params.ofn_stats << endl;
  ofstream out(params.ofn_stats.c_str());
  massert(out.is_open(), "could not open file %s", params.ofn_stats.c_str());
  out << "read_count\t" << stats.read_count << endl;
  out << "site_count\t" << stats.site_count << endl;
  out << "site_genomic_count\t" << stats.site_genomic_count << endl;
  out << "segment_count\t" << stats.segment_count << endl;
  out << "short_segment_count\t" << stats.short_segment_count << endl;
  out << "fend_count\t" << stats.fend_count << endl;
  out << "fend_genomic_count\t" << stats.fend_genomic_count << endl;
  out.close();

  return (0);
}
