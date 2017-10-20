#include <vector>
#include <string>
#include <string.h>
#include <map>

#include "util.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract Parser
////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum ParserType { ptString, ptFilename, ptInteger, ptDouble, ptBoolean };

struct Parser {
  // string identifier on command line
  string id;

  // desciption on usage
  string desc;

  // if true use only for usage
  bool dummy;

Parser(string _desc, bool _dummy) : desc(_desc), dummy(_dummy) {};
  // parse from string
  virtual void parse(char* arg) = 0;
  virtual ParserType type() = 0;

  // conversions
  virtual int to_int();
  virtual double to_double();
  virtual bool to_boolean();
  virtual string to_string();

  // class functions
  static string ParserType2String(ParserType ptype);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parsers
////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ParserInteger : public Parser {
  int value;
 ParserInteger(string _desc, int default_value=0, bool _dummy=false) : Parser(_desc, _dummy), value(default_value) {};
  ParserType type() { return ptInteger; };
  void parse(char* arg) { value = atoi(arg); };

  // conversions
  string to_string() { return std::to_string((long long int)value); };
  int to_int() { return value; };
};

struct ParserDouble : public Parser {
  double value;
 ParserDouble(string _desc, double default_value=0.0, bool _dummy=false) : Parser(_desc, _dummy), value(default_value) {};
  ParserType type() { return ptDouble; };
  void parse(char* arg) { value = atof(arg); };

  // conversions
  string to_string() { return std::to_string((long double)value); };
  double to_double() { return value; };
};

struct ParserBoolean : public Parser {
  bool value;
 ParserBoolean(string _desc, bool default_value=false, bool _dummy=false) : Parser(_desc, _dummy), value(default_value) {};
  ParserType type() { return ptBoolean; };
  void parse(char* arg) {
    massert(strlen(arg) == 1 && (arg[0] == 'T' || arg[0] == 'F'), "boolean must be T|F: %s", arg);
    value = (arg[0] == 'T');
  };

  // conversions
  string to_string() { return (value ? "T" : "F"); };
  bool to_boolean() { return value; };
};

struct ParserString : public Parser {
  string value;
 ParserString(string _desc, string default_value="", bool _dummy=false) : Parser(_desc, _dummy), value(default_value) {};
  ParserType type() { return ptString; };
  void parse(char* arg) { value = string(arg); };

  // conversions
  string to_string() { return value; };
};

struct ParserFilename : public ParserString {
 ParserFilename(string _desc, string default_value="", bool _dummy=false) : ParserString(_desc, default_value, _dummy) {};
  ParserType type() { return ptFilename; };
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Parameters {
 private:
  // raw values
  map<string, vector<char*> > m_values;

  // parsers
  map<string, Parser*> m_parsers;
  vector<Parser*> m_add_order;

  map<string, bool> m_mandatory;
  map<string, bool> m_set;

  // parse all into parameters
  Parser* get_parser(string id);
  char* get_value(string id);

 public:
  void usage(const char* name);

  // read from command line
  void read(int argc, char **argv);

  // parse values
  void add_parser(string id, Parser* param, bool mandatory=false);
  void parse(bool ignore_missing=false);

  // print all parsed values
  void print(ostream &os);

  // verify all mandatory params
  void verify_mandatory();

  const int get_int(string id) { return get_parser(id)->to_int(); };
  const double get_double(string id) { return get_parser(id)->to_double(); };
  const string get_string(string id) { return get_parser(id)->to_string(); };
  const bool get_bool(string id) { return get_parser(id)->to_boolean(); };
};
