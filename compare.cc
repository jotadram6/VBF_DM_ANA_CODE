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

  TFile* diff = new TFile("diff.root", "RECREATE");

  string filename1 = "test.root";
  string filename2 = "final.root";
  TFile* file1 = new TFile(filename1.c_str());
  TFile* file2 = new TFile(filename2.c_str());
  
  TDirectory* dir1 = (TDirectory*)file1->Get("NGenMuon");
  TDirectory* dir2 = (TDirectory*)file2->Get("DiJetCombinationsNmin");
  int count = 0;
  int all = 0;
  TIter iter(dir1->GetListOfKeys());
  TKey* key;
  while((key = (TKey*)iter())) {
    all++;
    string histname = key->GetName();
    
    if(dynamic_cast<TH2F*>(dir1->FindObjectAny(histname.c_str())) != NULL) {
      //      cout << histname << endl;
      TH2F* histo1 = (TH2F*)( dir1->FindObjectAny(histname.c_str()) );
      histname += "_0";
      if(switchNames(histname) != "no switch") histname = switchNames(histname);
      TH2F* histo2 = (TH2F*)( dir2->FindObjectAny(histname.c_str()) );
    
      if(histo1 != 0 && histo2 != 0) {
	bool passed = true;
	for(int i =0; i < (histo1->GetNbinsX()+1)*(histo1->GetNbinsX()+2); i++) {
	  if(abs(histo1->GetBinContent(i) - histo2->GetBinContent(i)) > 5) {
	    //TH1F* diffhist = new TH1F(histname.c_str(), histname.c_str(),  histo1->GetNbinsX(), histo1->GetXaxis()->GetXmin(), histo1->GetXaxis()->GetXmax());
	    // for(int j =0; j < histo1->GetNbinsX()+2; j++) {
	    //   diffhist->SetBinContent(j, histo1->GetBinContent(j) - histo2->GetBinContent(j));
	    // }
	    // diff->cd("");
	    // diffhist->Write();


	    cout << histname << " has different entries and is 2D" << endl;
	    count++;
	    passed = false;
	    break;
	  }
	}
      
      
      } else {
	cout << histname << " not found" << endl;
      }
      
    } else {
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
	for(int i =1; i < histo1->GetNbinsX()+1; i++) {
	  if(abs(histo1->GetBinContent(i) - histo2->GetBinContent(i)) > 500) {
	    TH1F* diffhist = new TH1F(histname.c_str(), histname.c_str(),  histo1->GetNbinsX(), histo1->GetXaxis()->GetXmin(), histo1->GetXaxis()->GetXmax());
	    for(int j =0; j < histo1->GetNbinsX()+2; j++) {
	      diffhist->SetBinContent(j, histo1->GetBinContent(j) - histo2->GetBinContent(j));
	    }
	    diff->cd("");
	    diffhist->Write();


	    cout << histname << " has different entries " << endl;
	    count++;
	    passed = false;
	    break;
	  }
	}
      
      
      } else {
	cout << histname << " not found" << endl;
      }
    }
  }
  cout << "Didn't match " << count << " out of " << all << " graphs" << endl;
  file1->Close();
  file2->Close();
  diff->Close();

  return count;
}


string switchNames(string name) {
  if(name == "LeadingJetMass_0") return "LeadingJetsMass_0";
  else if(name == "NTauJet1_0") return "NTau1_0";
  else if(name == "NTauJet2_0") return "NTau2_0";
  else if(name == "NVertices_0") return "NPVertices_0";

  else if(name == "Muon1MetDeltaPhiVsDiMuonCosDphi_0") return "Muon1MetDeltaPhiVsMuon1Muon2CosDphi_0";
  else if(name == "Tau1MetDeltaPhiVsDiTauCosDphi_0") return "Tau1MetDeltaPhiVsTau1Tau2CosDphi_0";
  else if(name == "DiMuon_Muon1DiJetDeltaPhi_0") return "Muon1Muon2_Muon1DiJetDeltaPhi_0";
  else if(name == "DiMuon_Muon2DiJetDeltaPhi_0") return "Muon1Muon2_Muon2DiJetDeltaPhi_0";
  else if(name == "DiMuonDeltaR_0") return "Muon1Muon2DeltaR_0";
  else if(name == "DiMuonDeltaPtDivSumPt_0") return "Muon1Muon2DeltaPtDivSumPt_0";
  else if(name == "DiMuonDeltaPt_0") return "Muon1Muon2DeltaPt_0";
  else if(name == "DiMuonOSLS_0") return "Muon1Muon2OSLS_0";
  else if(name == "DiMuonCosDphi_0") return "Muon1Muon2CosDphi_0";

  else if(name == "DiMuon_Muon1IsZdecay_0") return "Muon1Muon2_Muon1IsZmm_0";
  else if(name == "DiMuon_Muon2IsZdecay_0") return "Muon1Muon2_Muon2IsZmm_0";
  else if(name == "Muon1Tau1_Muon1IsZdecay_0") return "Muon1Tau1_Muon1IsZmm_0";
  else if(name == "Muon1Tau2_Muon1IsZdecay_0") return "Muon1Tau2_Muon1IsZmm_0";
  else if(name == "Muon2Tau1_Muon2IsZdecay_0") return "Muon2Tau1_Muon2IsZmm_0";
  else if(name == "Muon2Tau2_Muon2IsZdecay_0") return "Muon2Tau2_Muon2IsZmm_0";

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
