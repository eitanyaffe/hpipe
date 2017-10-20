#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <cmath>

#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <ctime>
#include <unordered_map>

using namespace std;

typedef uint64_t htype;

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

inline const htype hash_f(const string s)
{
  htype hash = 7;
  for (unsigned int i = 0; i < s.length(); i++)
    hash = hash*31 + s[i];
  return hash;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Parsing and handling user arguments
//////////////////////////////////////////////////////////////////////////////////////////////////

class UserParams {
public:
  vector<string> ifns1;
  vector<string> ifns2;
  string ofn1, ofn2;
  string mfn;
  string sfn;
};

void usage(const char* name, UserParams& params)
{
  fprintf(stderr, "usage: %s [options]\n", name);
  cout << " -ifn1 <fn>: input fasta read1 file, can have more than one" << endl;
  cout << " -ifn2 <fn>: input fasta read2 file" << endl;
  cout << " -ofn1 <fn>: output fasta read1" << endl;
  cout << " -ofn2 <fn>: output fasta read2" << endl;
  cout << " -mfn <fn>: multiplexcity output table" << endl;
  cout << " -sfn <fn>: stats output table" << endl;
  fprintf(stderr, "example: %s -ifn1 A1 -ifn2 A2 -ifn1 B1 -ifn2 B2 -ofn1 O1 -ofn2 O2 -mfn m -sfn s\n", name);
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

      if (option == "-ifn1")
	params.ifns1.push_back(arg);
      else if (option == "-ifn2")
	params.ifns2.push_back(arg);
      else if (option == "-ofn1")
	params.ofn1 = arg;
      else if (option == "-ofn2")
	params.ofn2 = arg;
      else if (option == "-mfn")
	params.mfn = arg;
      else if (option == "-sfn")
	params.sfn = arg;
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

inline string sort_concat_strings(const string& seq1, const string& seq2, const char delim)
{
  if (seq1.size() < seq2.size())
    return (seq1 + delim + seq2);
  if (seq1.size() > seq2.size())
    return (seq2 + delim + seq1);
  for (unsigned int i=0; i<seq1.size(); i++) {
    if (seq1[i] < seq2[i])
      return (seq1 + delim + seq2);
    if (seq1[i] > seq2[i])
      return (seq2 + delim + seq1);
  }
  return (seq1 + delim + seq2);
}

void traverse_fasta(ifstream& in1, ifstream& in2, ofstream& out1, ofstream& out2, unordered_map<htype, int>& multi,
		    double size, int& counter_total, int& counter_dups)
{
  const int MAXLINE = 1024;
  long counter = 0;
  char line1[4*MAXLINE], line2[4*MAXLINE];

  //  char tline1[MAXLINE], tline2[MAXLINE];

  clock_t begin = clock();
  clock_t time = begin;

  double total_seqs = 0;
  while(1) {
    int index = counter % 4;
    // cout << "line=" << counter+1 << " i=" << index << endl;

    in1.getline(line1 + index*MAXLINE, MAXLINE-1);
    in2.getline(line2 + index*MAXLINE, MAXLINE-1);

    bool eof1 = in1.eof();
    bool eof2 = in2.eof();
    massert(eof1 == eof2, "Error: eof of read1 and read2 files is not identical");
    bool eof = eof1;

    if (!eof && (!in1.good() || !in2.good())) {
      cerr << "Error reading line: " << counter+1
	   << ", rdstate1=" << in1.rdstate()
	   << ", rdstate2=" << in2.rdstate() << endl;
      if ( (in1.rdstate() & std::ifstream::failbit ) != 0  || (in2.rdstate() & std::ifstream::failbit ) != 0)
	cerr << "Probably line in file exceeds allocated maximal length: " << MAXLINE << endl;
      exit(-1);
    }

    //    if (strlen(tline1) > 0) {
    //  strcpy(line1 + index*MAXLINE, tline1);
    //  strcpy(line2 + index*MAXLINE, tline2);
    // }

    // check if sequence already encountered
    if (index == 3 && counter) {
      string seq1(line1 + MAXLINE);
      string seq2(line2 + MAXLINE);
      htype key = hash_f(sort_concat_strings(seq1, seq2, '_'));
      if (multi.find(key) == multi.end()) {
	multi[key] = 0;
	for (int i=0; i<4; i++) {
	  out1 << line1 + i*MAXLINE << endl;
	  out2 << line2 + i*MAXLINE << endl;
	}
      } else {
	counter_dups++;
      }
      multi[key]++;
      counter_total++;
    }

    if (eof)
      break;

    counter++;
    if (counter % 10000000 == 0) {
      clock_t ntime = clock();
      double secs = double(ntime - time) / CLOCKS_PER_SEC;
      time = ntime;
      cout << "progress: " << 100*((double)counter/4)/total_seqs << "%, delta sec=" << secs << endl;
    }

    // report once the estimated number of lines (given all sequences are same length)
    if (counter == 4) {
      ifstream::pos_type pos = in1.tellg();
      total_seqs = round(size/pos);
      cout << "estimated number of sequences (given length of first sequence): " << round(total_seqs/100000)/10 << "M" << endl;
    }
  }

  clock_t end = clock();
  double secs = double(end - begin) / CLOCKS_PER_SEC;
  cout << "total time traversing fasta file: " << secs/60 << " minutes" << endl;
}

void save_multi_table(unordered_map<htype, int>& multi, string ofn)
{
  // counter per each multiplicity
  map<int, int> table;

  cout << "computing multiplexity table..." << endl;
  for (unordered_map<htype, int>::iterator it = multi.begin(); it !=  multi.end(); ++it) {
    int count = (*it).second;
    if (table.find(count) == table.end())
      table[count] = 0;
    table[count]++;
  }

  cout << "saving multiplexity table to file: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  out << "multi" << "\t" << "count" << endl;

  for (map<int, int>::iterator it = table.begin(); it !=  table.end(); ++it) {
    int multi_value = (*it).first;
    int multi_count = (*it).second;
    out << multi_value << "\t" << multi_count << endl;
  }

  out.close();
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  UserParams params;
  parse_user_arguments(argc, argv, params);

  massert(params.ifns1.size() == params.ifns2.size(), "must have equal number of ifn1 and ifn2 input files");

  unordered_map<htype, int> multi;

  cout << "output read1: " << params.ofn1 << endl;
  cout << "output read2: " << params.ofn2 << endl;

  ofstream out1(params.ofn1.c_str());
  massert(out1.is_open(), "could not open file %s", params.ofn1.c_str());

  ofstream out2(params.ofn2.c_str());
  massert(out2.is_open(), "could not open file %s", params.ofn2.c_str());

  int counter_total=0, counter_dups=0;

  for (unsigned int i=0; i<params.ifns1.size(); i++) {
    string ifn1 = params.ifns1[i];
    string ifn2 = params.ifns2[i];

    cout << "input read1: " << ifn1 << endl;
    cout << "input read2: " << ifn2 << endl;

    double size1 = filesize(ifn1);

    ifstream in1(ifn1.c_str());
    massert(in1.is_open(), "could not open file %s", ifn1.c_str());

    ifstream in2(ifn2.c_str());
    massert(in2.is_open(), "could not open file %s", ifn2.c_str());
    traverse_fasta(in1, in2, out1, out2, multi, size1, counter_total, counter_dups);

    in1.close();
    in2.close();
  }
  cout << "total sequences: " << counter_total << endl;
  cout << "dup sequences: " << counter_dups << endl;
  cout << "yield (percentage of kept reads): " << (double) 100 * (counter_total-counter_dups) / counter_total << "%" << endl;

  out1.close();
  out2.close();

  save_multi_table(multi, params.mfn);

  // stats
  cout << "saving stats table to file: " << params.sfn << endl;
  ofstream out(params.sfn.c_str());
  massert(out.is_open(), "could not open file %s", params.sfn.c_str());
  out << "input" << "\t" << "output" << endl;
  out << counter_total << "\t" << counter_total-counter_dups << endl;
  out.close();

  return (0);
}
