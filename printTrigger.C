#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include "TBranch.h"
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
  if(argc < 2) {
    cout << "No file provided" << endl;
    exit(1);
  } 
  string file = argv[1];
  if(file.find(".root") == string::npos) {
    cout << "Input file is not a root file" << endl;
    exit(1);
  } else if(argc > 2) {
    cout << "More than one file provided, these files will be ignored:" << endl;
    for(int i=2; i < argc; i++) {
      cout << argv[i] << endl;
    } 
    cout << endl;
  }

  TFile* f = TFile::Open(file.c_str());
  f->cd("TNT");
  TTree* BOOM = (TTree*)f->Get("TNT/BOOM");

  vector<string>* Trigger_names = 0;
  BOOM->SetBranchAddress("Trigger_names", &Trigger_names);
  BOOM->GetEntry(0);

  cout << "Triggers: " << endl << endl;
  for(int i=0; i < (int)Trigger_names->size(); i++) {
    cout << Trigger_names->at(i) << endl;
  }
  f->Close();

}


