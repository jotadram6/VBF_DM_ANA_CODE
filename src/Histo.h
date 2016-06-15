#ifndef Histo_h
#define Histo_h

// system include files
#include <memory>

// user include files
#include <Math/VectorUtil.h>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
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
#include "Cut_enum.h"

using namespace std;

class Histogramer {

 public:
  Histogramer(int, std::string, std::string, std::string);
  ~Histogramer();
  int NFolders;
  int Npdf;
  
  vector<int> get_folders();
  int get_Nhists();
  double get_start(int);
  double get_width(int);
  int get_nbins(int);
  map<string,pair<int,int>>* get_cuts();
  
 private:
  TFile * outfile;
  map<string, pair<int,int>> cuts;
  vector<string> folders;

  std::map<string, std::vector<TH1*> > Generator_Histogram;
  vector<pair<string, std::array<double, 3> > > Generator_info;

  int index(int, int);
  void write_histogram();
  void read_hist(string);
  void read_cuts(string);
  void fill_histogram(TObject*);

  std::unordered_map<string, CUTS> cut_num = { {"NGenTau", CUTS::eGTau}, {"NGenTop", CUTS::eGTop}, {"NGenElectron", CUTS::eGElec}, \
					       {"NGenMuon", CUTS::eGMuon}, {"NGenZ", CUTS::eGZ}, {"NGenW", CUTS::eGW}, {"NGenHiggs", CUTS::eGHiggs}, \
					       {"NRecoVertex", CUTS::eRVertex}, {"NRecoMuon1", CUTS::eRMuon1}, {"NRecoMuon2", CUTS::eRMuon2}, \
					       {"NRecoElectron1", CUTS::eRElec1}, {"NRecoElectron2",CUTS::eRElec2}, {"NRecoTau1", CUTS::eRTau1},  \
					       {"NRecoTau2", CUTS::eRTau2}, {"NRecoJet1", CUTS::eRJet1}, {"NRecoJet2", CUTS::eRJet2}, \
					       {"NRecoCentralJet", CUTS::eRCenJet}, {"NRecoBJet", CUTS::eRBJet}, {"NRecoTriggers1", CUTS::eRTrig1},
					       {"NRecoTriggers2", CUTS::eRTrig2}, {"NRecoFirstLeadingJet", CUTS::eR1stJet}, {"NRecoSecondLeadingJet", CUTS::eR2ndJet},
					       {"NRecoMuon1MetTopology", CUTS::eTMuon1}, {"NRecoMuon2MetTopology", CUTS::eTMuon2},
					       {"NRecoElectron1MetTopology", CUTS::eTElec1}, {"NRecoElectron2MetTopology", CUTS::eTElec2},
					       {"NRecoTau1MetTopology", CUTS::eTTau1}, {"NRecoTau2MetTopology", CUTS::eTTau2}, {"NDiMuonCombinations", CUTS::eDiMuon},
					       {"NDiElectronCombinations", CUTS::eDiElec}, {"NDiTauCombinations", CUTS::eDiTau}, {"NDiJetCombinations", CUTS::eDiJet},
					       {"NMuon1Tau1Combinations", CUTS::eMuon1Tau1}, {"NMuon1Tau2Combinations", CUTS::eMuon1Tau2},
					       {"NMuon2Tau1Combinations", CUTS::eMuon2Tau1}, {"NMuon2Tau2Combinations", CUTS::eMuon2Tau2},
					       {"NElectron1Tau1Combinations", CUTS::eElec1Tau1}, {"NElectron1Tau2Combinations", CUTS::eElec1Tau2},
					       {"NElectron2Tau1Combinations", CUTS::eElec2Tau1}, {"NElectron2Tau2Combinations", CUTS::eElec2Tau2},
					       {"NSusyCombinations", CUTS::eSusyCom} };

};

#endif

