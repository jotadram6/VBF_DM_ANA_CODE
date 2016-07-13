#include "Analyzer.h"

int main (int argc, char* argv[]) {
  
  if(argc < 3) {
    std::cout << "You have entered too little arguments, please type:" << std::endl;
    std::cout << "./Analyzer infile.root outfile.root" << std::endl;
    exit(EXIT_FAILURE);
    ////TODO, add ability to give multiple input files
  } else if(argc > 3) {
    std::cout << "You have entered too many arguments, please type:" << std::endl;
    std::cout << "./Analyzer infile.root outfile.root" << std::endl;
    exit(EXIT_FAILURE);
  }

  string name = argv[1];
  ifstream ifile(name);
  if ( !ifile && name.find("root://") == string::npos && name.find("root\\://") == string::npos) {
    std::cout << "The file '" << argv[1] << "' doesn't exist" << std::endl;
    exit(EXIT_FAILURE);
  }



  Analyzer testing(argv[1], argv[2]);
  //  Analyzer testing("TNT.root", "test.root");

  for(int i=0; i < 10000; i++) {
    testing.clear_values();
    testing.preprocess(i);
    testing.fill_histogram();
  }
  testing.printCuts();
  return 0;
}
