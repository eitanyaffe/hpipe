#include "Params.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
  Parameters params;

  // add parsers
  params.add_parser("F", new ParserFilename("filename", "x.cpp"), true);
  params.add_parser("S", new ParserString("name", "ABC"));
  params.add_parser("D", new ParserDouble("fraction", 0.5));
  params.add_parser("B", new ParserBoolean("flag", true));
  params.add_parser("N", new ParserInteger("size", 0));
  params.add_parser("Vi", new ParserInteger("field i, i=1...N", 0));

  if (argc == 1) {
    params.usage(argv[0]);
    exit(0);
  }

  // read command line params
  params.read(argc, argv);
  params.parse(true);

  params.verify_mandatory();

  int N = params.get_int("N");
  for (int i=1; i<=N; i++) {
    string s = to_string((long long)i);
    params.add_parser(string("V")+s, new ParserInteger(string("field ")+s, 0));
  }

  params.parse(false);
  params.print(cerr);

  string F = params.get_string("F");
  string S = params.get_string("S");
  double D = params.get_double("D");
  double B = params.get_bool("B");
  cout << "F=" << F << " S=" << S << " D=" << D << " N=" << N << " B=" << B << endl;
}
