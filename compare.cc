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

string switchNames(string);

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
    if(switchNames(histname) != "no switch") histname = switchNames(histname);
    TH1F* histo2 = (TH1F*)( dir2->FindObjectAny(histname.c_str()) );
    
    if(histo1 != 0 && histo2 != 0) {
      bool passed = true;
      if(histo1->GetNbinsX() != histo2->GetNbinsX()) {
	cout << histname << " has " << histo1->GetNbinsX() << " and " << histo2->GetNbinsX() << " bins" << endl;
	count++;
	continue;
      }
      for(int i =0; i < histo1->GetNbinsX(); i++) {
	if(abs(histo1->GetBinContent(i) - histo2->GetBinContent(i)) > 1) {

	    cout << histname << " has different entries " << endl;
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


string switchNames(string name) {
  if(name == "LeadingJetMass_0") return "LeadingJetsMass_0";
  else if(name == "LeadingJetPt_0") return "LeadingJetsPt_0";
  else if(name == "LeadingJetDeltaR_0") return "LeadingJetsDeltaR_0";
  else if(name == "LeadingJetDeltaEta_0") return "LeadingJetsDeltaEta_0";
  else if(name == "DiTau_Tau1DiJetDeltaPhi_0") return "Tau1Tau2_Tau1DiJetDeltaPhi_0";
  else if(name == "DiTau_Tau2DiJetDeltaPhi_0") return "Tau1Tau2_Tau2DiJetDeltaPhi_0";
  else if(name == "DiTauDeltaR_0") return "Tau1Tau2DeltaR_0";
  else if(name == "DiTauDeltaPtDivSumPt_0") return "Tau1Tau2DeltaPtDivSumPt_0";
  else if(name == "DiTauDeltaPt_0") return "Tau1Tau2DeltaPt_0";
  else if(name == "DiTauOSLS_0") return "Tau1Tau2OSLS_0";
  else if(name == "DiTauCosDphi_0") return "Tau1Tau2CosDphi_0";
  else return "no switch";
}
