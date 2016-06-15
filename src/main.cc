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


  ifstream ifile(argv[1]);
  if ( !ifile) {
    std::cout << "The file '" << argv[1] << "' doesn't exist" << std::endl;
    exit(EXIT_FAILURE);
  }



  Analyzer* testing = new Analyzer(argv[1], argv[2]);

  for(int i=0; i < testing->nentries; i++) {
    testing->clear_values();
    testing->preprocess(i);
  }
  testing->printCuts();
}
