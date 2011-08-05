/**** written by Thomas Schoenemann as a private person without employment, September 2009 ***/

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

typedef unsigned int uint;

int main(int argc, char** argv) {

  if (argc == 1) {

    std::cout << "ngrep <string> [<filename 1>] [<filename 2>] ... " << std::endl << std::endl 
	      << "     a program that filters lines from text-files: " << std::endl
	      << "     it prints all lines that do not contain the first argument." << std::endl
	      << "     Hence, it is the complement to the well-known \"grep\"-command" << std::endl;

    exit(0);
  }
    
  std::string key = argv[1];

  std::vector<std::string> filenames;
  if (argc == 2)
    filenames.push_back("&0");
  else {
    for (int k=2; k < argc; k++)
      filenames.push_back(std::string(argv[k]));
  }

  char cline[65536];
  for (uint i=0; i < filenames.size(); i++) {

    std::ifstream finstream;
    if (filenames[i] != "&0") {
      finstream.open(filenames[i].c_str());

      if (!finstream.is_open()) {
	std::cerr << "ERROR: could not open file \"" << filenames[i] << "\". Exiting" << std::endl;
	exit(1);
      }
    }

    std::istream& instream = (filenames[i] != "&0") ? finstream : std::cin;

    while (instream.getline(cline,65536)) {

      std::string line = cline;
      if (line.find(key) > line.size())
	std::cout << line << std::endl;
    }
    
    if (filenames[i] != "&0")
      finstream.close();
  }

}
