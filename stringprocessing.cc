/*** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include "stringprocessing.hh"
#include "makros.hh"

bool is_whitespace(char c) {

  if (c == ' ' || c == '\n' || c == '\t' || c == 10 || c == '\r')
    return true;
  return false;
}

bool is_uppercase(char c) {

  if (c >= 'A' && c <= 'Z')
    return true;
//   if (c == 'Ä' || c == 'Ö' || c == 'Ü')
//      return true;
  
  return false;
}

char downcase(char c) {

  if (c >= 'A' && c <= 'Z')
    return c + ('a' - 'A');

  //outcommented since gcc claims these are several characters
//   if (c == 'Ä')
//     return 'ä';
//   if (c == 'Ö')
//     return 'ö';
//   if (c == "Ü")
//     return "ü";
//   if (c == 'É')
//     return 'é'

  return c;
}

std::string downcase(std::string s) {

  std::string ls=s;
  for (uint k=0; k < ls.size(); k++)
    ls[k] = downcase(ls[k]);

  return ls;
}

bool is_natural_number(const std::string s) {

  for (uint i=0; i < s.size(); i++) {
    if (s[i] < '0' || s[i] > '9')
      return false;
  }
  return true;
}


//to avoid the awkward strcmp routine
bool strings_equal(std::string s1, std::string s2) {
    return (s1 == s2);
}

void tokenize(const std::string& s, std::vector<std::string>& tokens, char separator, bool empty_tokens) {

  bool last_sep = true;
  tokens.clear();

  for (uint i=0; i < s.size(); i++) {

    if (s[i] == separator) {
      
      if (empty_tokens)
	tokens.push_back(std::string());

      last_sep = true;
    }
    else {
      if (last_sep)
	tokens.push_back(std::string());
      tokens.back() += s[i];
      last_sep = false;
    }
  }
}

void tokenize_with_stringsep(const std::string& org_s, std::vector<std::string>& tokens, std::string sep_string, bool empty_tokens) {

  tokens.clear();

  std::string s = org_s;

  while (true) {

    int pos = s.find(sep_string);

    if (pos > s.size()) {
      if (s != "")
	tokens.push_back(s);
      break;
    }
    else {
      if (empty_tokens || pos != 0)
	tokens.push_back(s.substr(0,pos));

      s = s.substr(pos+sep_string.size());
    }
  }
}

bool string_ends_with(std::string s, std::string suffix) {

  if (s.size() < suffix.size())
    return false;
  else
    return s.substr(s.size() - suffix.size(), suffix.size()) == suffix;
}
