#ifndef Analyzer_h
#define Analyzer_h

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
#include "Particle.h"
#include "Histo.h"


using namespace std;

class Analyzer {

 public:
  Analyzer(string, string);
  void clear_values();
  void preprocess(int);
  int fillCuts(map<string,pair<int,int>>*);
  void printCuts();
  int nentries;

  
 private:
  void getInputs();
  void setupJob(string);
  void initializePileupInfo(string, string);  
  void read_info(string);
  void setupGeneral(TTree*, string);

  void SmearLepton(Lepton*, CUTS, PartStats*);
  void SmearJet();

  bool JetMatchesLepton(Lepton*, TLorentzVector, double, CUTS);
  TLorentzVector* matchLeptonToGen(TLorentzVector*, Lepton*, CUTS);
  TLorentzVector* matchTauToGen(TLorentzVector*, double);

  void ExtractNumberOfTauNu();  
  void ExtractNumberOfGoodGen(int, int, CUTS, PartStats*);
  void ExtractNumberOfGoodReco(Lepton*, CUTS, CUTS, PartStats*);
  void ExtractNumberOfGoodRecoJets(CUTS, PartStats*);

  void passRecoLeptonMetTopologyCuts(Lepton*, CUTS,CUTS, PartStats*);
  void passLeptonComboTopologyCut(Lepton*, Lepton*, CUTS,CUTS,CUTS, PartStats*);
  void passDiJetTopologyCuts(PartStats*);

  void SusyTopologyCuts();
  bool passTriggerCuts(string);

  double CalculateLeptonMetMt(TLorentzVector);
  double DiParticleMass(TLorentzVector, TLorentzVector, string);
  bool passDiParticleApprox(TLorentzVector, TLorentzVector, string);
  bool isZdecay(TLorentzVector, Lepton*);

  bool isOverlaping(TLorentzVector, Lepton*, CUTS, double);
  bool passProng(string, int);
  bool isInTheCracks(float);
  bool passedLooseJetID(int);

  double CalculatePZeta(TLorentzVector, TLorentzVector);
  double CalculatePZetaVis(TLorentzVector, TLorentzVector);
  double normPhi(double);
  double absnormPhi(double);

  void updateMet();
  double getPileupWeight(float);


  Generated* _Gen;
  Electron* _Electron;
  Muon* _Muon;
  Taus* _Tau;
  Jet* _Jet;
  Histogramer* histo;

  unordered_map<string, PartStats> distats;
  unordered_map<string, pair<int,int> > prevTrig;
  std::array<std::vector<int>, static_cast<int>(CUTS::enumSize)> goodParts;
  
  TFile* f;
  TTree* BOOM;

  vector<int> cuts_per, cuts_cumul;

  TLorentzVector theMETVector;
  double deltaMEx, deltaMEy, sumpxForMht, sumpyForMht, sumptForHt, phiForMht;
  double againstElectron, againstMuon, maxIso, minIso;
  int leadIndex;
  bool isData, CalculatePUSystematics;

  vector<double>* Trigger_decision = 0;
  vector<string>* Trigger_names = 0;
  float nTruePU = 0;
  int bestVertices = 0;
  double Met_px = 0;
  double Met_py = 0;
  double Met_pz = 0;
  
  TH1F *hPUmc = new TH1F("hPUmc", "hPUmc", 100, 0, 100);
  TH1F *hPUdata = new TH1F("hPUdata", "hPUdata", 100, 0, 100);
  double pu_weight;
  
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
