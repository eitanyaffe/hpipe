#ifndef __UTIL__
#define __UTIL__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

enum InputType { itFastq, itFasta };

void split_string(string str, char delim, vector<string>& result);
void split_line(istream &in, vector<string> &fields, char delim);
void massert(bool cond, char *fmt, ...);
void mexit(char *fmt, ...);
int get_field_index(string field, const vector<string>& titles);

string reverse_complement(string str);

#endif // __UTIL__
