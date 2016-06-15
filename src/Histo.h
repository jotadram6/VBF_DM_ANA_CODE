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

class Histogramer {

 public:
  Histogramer(int, string);
  ~Histogramer();
  int NFolders;
  int Npdf;
  
  vector<string> get_cuts();
  vector<int> get_folders();
  int get_Nhists();
  double get_start(int);
  double get_width(int);
  int get_nbins(int);
  
 private:
  TFile * outfile;
  vector<string> cuts;
  vector<string> folders;
  vector<pair<int,int>> cut_range;

  std::map<string, std::vector<TH1*> > Generator_Histogram;
  vector<pair<string, std::array<double, 3> > > Generator_info;

  int index(int, int);
  void write_histogram();
  void read_hist();
  void read_cuts();
  void fill_histogram(TObject*);

  std::map<string, CUTS> Cut_num = { {"NGenTau", CUTS::eGTau}, {"NGenTop", CUTS::eGTop}, {"NGenElectron", CUTS::eGElec}, \
    {"NGenMuon", CUTS::eGMuon}, {"NGenZ", CUTS::eGZ}, {"NGenW", CUTS::eGW}, {"NGenHiggs", CUTS::eGHiggs}, \
    {"NREcoVertex", CUTS::eRVertex}, {"NRecoMuon1", CUTS::eRMuon1}, {"NRecoMuon2", CUTS::eRMuon2}, \
    {"NRecoElectron1", CUTS::eRElec1}, {"NRecoElectron2",CUTS::eRElec2}, {"NRecoTau1", CUTS::eRTau1}, \
				     {"NRecoTau2", CUTS::eRTau2}, {"NRecoJet1", CUTS::eRJet1}, {"NRecoJet2", CUTS::eRJet2}, {"NRecoCentralJet", CUTS::eRCenJet} };

};
#endif
UTS::eRCenJet}, 
				     {"NRecoBJet", CUTS::eRBJet }
{"NRecoTriggers1", CUTS::
{"NRecoTriggers2", CUTS::
{"NRecoFirstLeadingJet", CUTS::
{"NRecoSecondLeadingJet", CUTS::

{"NRecoMuon1MetTopology", CUTS::
{"NRecoMuon2MetTopology", CUTS::
{"NRecoElectron1MetTopology", CUTS::
{"NRecoElectron2MetTopology", CUTS::
{"NRecoTau1MetTopology", CUTS::
{"NRecoTau2MetTopology", CUTS::

{"NDiMuonCombinations", CUTS::
{"NDiElectronCombinations", CUTS::
{"NDiTauCombinations", CUTS::
{"NDiJetCombinations", CUTS::

{"NMuon1Tau1Combinations", CUTS::
{"NMuon1Tau2Combinations", CUTS::
{"NMuon2Tau1Combinations", CUTS::
{"NMuon2Tau2Combinations", CUTS::
{"NElectron1Tau1Combinations", CUTS::
{"NElectron1Tau2Combinations", CUTS::
{"NElectron2Tau1Combinations", CUTS::
{"NElectron2Tau2Combinations", CUTS::

{"NSusyCombinations", CUTS::


};
#endif
