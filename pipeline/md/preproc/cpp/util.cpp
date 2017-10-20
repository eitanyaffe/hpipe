#include "util.h"

string reverse_complement(string str)
{
  string result = str;
  int N = str.length();
  for (int i=0; i<N; i++) {
    char c = str[N-i-1], r;
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

void split_string(string str, char delim, vector<string>& result)
{
  unsigned int prev = 0;
  for (unsigned int i=0; i<str.size(); i++) {
    char c = str[i];
    if(c == delim) {
      result.push_back(str.substr(prev, i-prev));
      prev = i;
    }
  }
  if (prev < str.size())
    result.push_back(str.substr(prev, str.size()-prev));
}

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

void mexit(char *fmt, ...)
{
  fprintf(stderr, "Error: ");

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  fprintf(stderr, "\n");
  exit(-1);
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

