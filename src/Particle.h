#ifndef Particle_h
#define Particle_h

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
#include "tokenizer.hpp"
#include "Cut_enum.h"

using namespace std;

struct PartStats {
  unordered_map<string,double> dmap;
  unordered_map<string,string> smap;
  unordered_map<string,pair<double,double> > pmap;
  unordered_map<string,bool> bmap;
};


class Particle {

 public:
  Particle();
  Particle(TTree*, string, string);
  Particle& operator =(Particle&);
  virtual ~Particle() {};
  void getPartStats(string);

  vector<double>* pt = 0;
  vector<double>* eta = 0;
  vector<double>* phi = 0;
  vector<double>* energy = 0;

  unordered_map<string, PartStats> pstats;

  vector<TLorentzVector> smearP;
};

///template <typename T> vector<TLorentzVector> Particle<T>::smearP = {}; 


/////////////////////////////////////////////////////////////////
class Generated : public Particle {

public:
  Generated();
  Generated(TTree*, string);
  Generated& operator =(Generated&);

  vector<double>  *pdg_id = 0;
  vector<double>  *motherpdg_id = 0;
  vector<double>  *status = 0;
  vector<int>  *BmotherIndex = 0;

};

/////////////////////////////////////////////////////////////////////////
class Jet : public Particle {

public:
  Jet(TTree*, string);

  vector<double>* neutralHadEnergyFraction = 0;
  vector<double>* neutralEmEmEnergyFraction = 0;
  vector<int>*    numberOfConstituents = 0;
  vector<double>* muonEnergyFraction = 0;
  vector<double>* chargedHadronEnergyFraction = 0;
  vector<int>*    chargedMultiplicity = 0;
  vector<double>* chargedEmEnergyFraction = 0;
  vector<int>*    partonFlavour = 0;
  vector<double>* bDiscriminator = 0;
};


class Lepton : public Particle {

public: 
  Lepton(TTree*, string, string);

  vector<double>* charge = 0;

};

class Electron : public Lepton {

public:
  Electron(TTree*, string);

   vector<int>     *isPassVeto = 0;
   vector<int>     *isPassLoose = 0;
   vector<int>     *isPassMedium = 0;
   vector<int>     *isPassTight = 0;
   vector<int>     *isPassHEEPId = 0;
   vector<double>  *isoChargedHadrons = 0;
   vector<double>  *isoNeutralHadrons = 0;
   vector<double>  *isoPhotons = 0;
   vector<double>  *isoPU = 0;
};

class Muon : public Lepton {

public: 
  Muon(TTree*, string);

   vector<bool>* tight = 0;
   vector<bool>* soft = 0;
   vector<double>* isoCharged = 0;
   vector<double>* isoNeutralHadron = 0;
   vector<double>* isoPhoton = 0;
   vector<double>* isoPU = 0;
};

class Taus : public Lepton {

 public:
  Taus(TTree*, string);

   vector<int>     *decayModeFindingNewDMs = 0;
   vector<double>  *nProngs = 0;
   pair<vector<int>*,vector<int>* > againstElectron = make_pair(nullptr,nullptr);
   pair<vector<int>*,vector<int>* > againstMuon = make_pair(nullptr,nullptr);
   pair<vector<int>*,vector<int>* > minIso = make_pair(nullptr,nullptr);
   pair<vector<int>*,vector<int>* > maxIso = make_pair(nullptr,nullptr);
   vector<double>  *leadChargedCandPt = 0;
};



#endif
