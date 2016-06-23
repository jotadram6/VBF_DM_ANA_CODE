#include <Math/VectorUtil.h>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCollection.h>
#include <TKey.h>
#include <TTree.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <utility>
#include <TROOT.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TLorentzVector.h>
#include <TEnv.h>

using namespace std;

int compare() {

  string filename1 = "test.root";
  string filename2 = "final.root";
  TFile* file1 = new TFile(filename1.c_str());
  TFile* file2 = new TFile(filename2.c_str());
  
  TDirectory* dir1 = (TDirectory*)file1->Get("NGenMuon");
  TDirectory* dir2 = (TDirectory*)file2->Get("DiJetCombinationsNmin");
  int count = 0;

  TIter iter(dir1->GetListOfKeys());
  TKey* key;
  while((key = (TKey*)iter())) {
    string histname = key->GetName();
    TH1F* histo1 = (TH1F*)( dir1->FindObjectAny(histname.c_str()) );
    histname += "_0";
    TH1F* histo2 = (TH1F*)( dir2->FindObjectAny(histname.c_str()) );
    
    if(histo1 != 0 && histo2 != 0) {
      bool passed = true;
      if(histo1->GetNbinsX() != histo2->GetNbinsX()) {
	cout << histname << " has " << histo1->GetNbinsX() << " and " << histo2->GetNbinsX() << " bins" << endl;
	count++;
	continue;
      }
      for(int i =0; i < histo1->GetNbinsX(); i++) {
	if(histo1->GetBinContent(i) != histo2->GetBinContent(i)) {
	  cout << histname << " has different entries" << endl;
	  count++;
	  passed = false;
	  break;
	}
      }
      if(passed) cout << "******" << histname << " passed**********" << endl;
    } else {
      cout << histname << " not found" << endl;
      count++;
    }
  }
  return count;
}
