#include "Analyzer.h"
#define ival(x) static_cast<int>(x)
/*For speed, uncomment these commands and the define statement in the h file.  It will lead to 
  unchecked bounds on values, but will improve speed */
//#define const
//#define at(x) operator[](x)
typedef vector<int>::iterator vec_iter;

//Filespace that has all of the .in files
const string FILESPACE = "PartDet/";

//////////PUBLIC FUNCTIONS////////////////////

///Constructor
Analyzer::Analyzer(string infile, string outfile) : hPUmc(new TH1F("hPUmc", "hPUmc", 100, 0, 100)), hPUdata(new TH1F("hPUdata", "hPUdata", 100, 0, 100)), histo(1, FILESPACE+"Hist_entries.in", FILESPACE+"Cuts.in", outfile) {

  f = TFile::Open(infile.c_str());
  f->cd("TNT");
  BOOM = (TTree*)f->Get("TNT/BOOM");
  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "setup" << std::endl;
  std::cout << nentries << std::endl;

  setupGeneral(BOOM,infile);

  isData = distats["Run"].bmap.at("isData");
  CalculatePUSystematics = distats["Run"].bmap.at("CalculatePUSystematics");

  prevTrig["Trigger1"] = make_pair(0,0);
  prevTrig["Trigger2"] = make_pair(0,0);
  initializePileupInfo(distats["Run"].smap.at("MCHistos"), distats["Run"].smap.at("DataHistos"));

  //////need to initialize histo and get values for cut arrays

  cuts_per.resize(histo.get_cuts()->size());
  cuts_cumul.resize(histo.get_cuts()->size());

  
  _Gen = new Generated(BOOM, FILESPACE + "Gen_info.in");
  _Electron = new Electron(BOOM, FILESPACE + "Electron_info.in");
  _Muon = new Muon(BOOM, FILESPACE + "Muon_info.in");
  _Tau = new Taus(BOOM, FILESPACE + "Tau_info.in");
  _Jet = new Jet(BOOM, FILESPACE + "Jet_info.in");
}

////destructor
Analyzer::~Analyzer() {
  delete f;
  delete _Gen;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
}


///resets values so analysis can start
void Analyzer::clear_values() {
  for(int i=0; i < (int)goodParts.size(); i++) {
    goodParts[i].clear();
  }
  deltaMEx=0;
  deltaMEy=0;
  sumpxForMht = 0.0;
  sumpyForMht = 0.0;
  sumptForHt  = 0.0;
  leadIndex=-1;
}

///Function that does most of the work.  Calculates the number of each particle
void Analyzer::preprocess(int event) {
  BOOM->GetEntry(event);

  //TODO: add in pdf vector(set to 1 for now);
  
  theMETVector.SetPxPyPzE(Met_px, Met_py, Met_pz, sqrt(pow(Met_px,2) + pow(Met_py,2)));
  pu_weight = (!isData && CalculatePUSystematics) ? getPileupWeight(nTruePU) : 1.0;

   // clock_t t1;
   // t1 = clock();
 
  // SET NUMBER OF GEN PARTICLES
  // TODOGeneralize to remove magic numbers
  if(!isData){
    PartStats genStat = _Gen->pstats["Gen"];
    getGoodGen(15, 2, CUTS::eGTau, genStat);
    getGoodGen(6, 2, CUTS::eGTop, genStat);
    getGoodGen(11, 1, CUTS::eGElec, genStat);
    getGoodGen(13, 1, CUTS::eGMuon, genStat);
    getGoodGen(23, 2, CUTS::eGZ, genStat);
    getGoodGen(24, 2, CUTS::eGW, genStat);
    getGoodGen(25, 2, CUTS::eGHiggs, genStat);
    getGoodTauNu();
  }

   // gent+= clock() - t1;
   // t1 = clock();

  //////Smearing  
  smearLepton(*_Electron, CUTS::eGElec, _Electron->pstats["Smear"]);
  smearLepton(*_Muon, CUTS::eGMuon, _Muon->pstats["Smear"]);
  smearLepton(*_Tau, CUTS::eGTau, _Tau->pstats["Smear"]);
  smearJet(_Jet->pstats["Smear"]);

  // smeart += clock() - t1;
  // t1 = clock();

  //////Triggers and Vertices
  goodParts[ival(CUTS::eRVertex)].resize(bestVertices);
  if(passTriggerCuts("Trigger1")) goodParts[ival(CUTS::eRTrig1)].resize(1);
  if(passTriggerCuts("Trigger2")) goodParts[ival(CUTS::eRTrig2)].resize(1);



  // trigt += clock() - t1;
  // t1 = clock(); 

  // // SET NUMBER OF RECO PARTICLES
  // // MUST BE IN ORDER: Muon/Electron, Tau, Jet
  getGoodRecoLeptons(*_Electron, CUTS::eRElec1, CUTS::eGElec, _Electron->pstats["Elec1"]);
  getGoodRecoLeptons(*_Electron, CUTS::eRElec2, CUTS::eGElec, _Electron->pstats["Elec2"]);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon1, CUTS::eGMuon, _Muon->pstats["Muon1"]);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon2, CUTS::eGMuon, _Muon->pstats["Muon2"]);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau1, CUTS::eGTau, _Tau->pstats["Tau1"]);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau2, CUTS::eGTau, _Tau->pstats["Tau2"]);

  // recolt += clock() - t1;
  // t1 = clock();

  getGoodRecoJets(CUTS::eRJet1, _Jet->pstats["Jet1"]);
  // jet1t += clock() - t1;
  // t1 = clock();

  getGoodRecoJets(CUTS::eRJet2, _Jet->pstats["Jet2"]);
  // jet2t += clock() - t1;
  // t1 = clock();

  getGoodRecoJets(CUTS::eRCenJet, _Jet->pstats["CentralJet"]);
  // cent += clock() - t1;
  // t1 = clock();

  getGoodRecoJets(CUTS::eRBJet, _Jet->pstats["BJet"]);
  // bjett += clock() - t1;
  // t1 = clock();

  getGoodRecoJets(CUTS::eR1stJet, _Jet->pstats["FirstLeadingJet"]);
  leadIndex = goodParts[ival(CUTS::eR1stJet)].at(0); 
  getGoodRecoJets(CUTS::eR2ndJet, _Jet->pstats["SecondLeadingJet"]);

  // leadt += clock() - t1;
  // t1 = clock();

  updateMet();

  /////  SET NUMBER OF RECO MET TOPOLOGY PARTICLES
  getGoodMetTopologyLepton(*_Electron, CUTS::eRElec1, CUTS::eTElec1, _Electron->pstats["Elec1"]);
  getGoodMetTopologyLepton(*_Electron, CUTS::eRElec2, CUTS::eTElec2, _Electron->pstats["Elec2"]);
  getGoodMetTopologyLepton(*_Muon, CUTS::eRMuon1, CUTS::eTMuon1, _Muon->pstats["Muon1"]);
  getGoodMetTopologyLepton(*_Muon, CUTS::eRMuon2, CUTS::eTMuon2, _Muon->pstats["Muon2"]);
  getGoodMetTopologyLepton(*_Tau, CUTS::eRTau1, CUTS::eTTau1, _Tau->pstats["Tau1"]);
  getGoodMetTopologyLepton(*_Tau, CUTS::eRTau2, CUTS::eTTau2, _Tau->pstats["Tau2"]);

 
  // mett += clock() - t1;
  // t1 = clock();

  ///VBF Susy cut on leadin jets
  if(goodParts[ival(CUTS::eR1stJet)].at(0) != -1 && goodParts[ival(CUTS::eR2ndJet)].at(0) != -1) VBFTopologyCut();

  // susyt += clock() - t1;
  // t1 = clock();

  /////lepton lepton topology cuts
  getGoodLeptonCombos(*_Tau, *_Electron, CUTS::eRTau1, CUTS::eRElec1, CUTS::eElec1Tau1, distats["Electron1Tau1"]);
  getGoodLeptonCombos(*_Tau, *_Electron, CUTS::eRTau1, CUTS::eRElec2, CUTS::eElec2Tau1, distats["Electron2Tau1"]);
  getGoodLeptonCombos(*_Tau, *_Muon, CUTS::eRTau1, CUTS::eRMuon1, CUTS::eMuon1Tau1, distats["Muon1Tau1"]);
  getGoodLeptonCombos(*_Tau, *_Muon, CUTS::eRTau1, CUTS::eRMuon2, CUTS::eMuon2Tau1, distats["Muon2Tau1"]);
  getGoodLeptonCombos(*_Tau, *_Electron, CUTS::eRTau2, CUTS::eRElec1, CUTS::eElec1Tau2, distats["Electron1Tau2"]);
  getGoodLeptonCombos(*_Tau, *_Electron, CUTS::eRTau2, CUTS::eRElec2, CUTS::eElec2Tau2, distats["Electron2Tau2"]);
  getGoodLeptonCombos(*_Tau, *_Muon, CUTS::eRTau2, CUTS::eRMuon1, CUTS::eMuon1Tau2, distats["Muon1Tau2"]);
  getGoodLeptonCombos(*_Tau, *_Muon, CUTS::eRTau2, CUTS::eRMuon2, CUTS::eMuon2Tau2, distats["Muon2Tau2"]);

  // combot += clock() - t1;
  // t1 = clock();

  ////DIlepton topology cuts
  getGoodLeptonCombos(*_Tau, *_Tau, CUTS::eRTau1, CUTS::eRTau2, CUTS::eDiTau, distats["DiTau"]);
  getGoodLeptonCombos(*_Electron, *_Electron, CUTS::eRElec1, CUTS::eRElec2, CUTS::eDiElec, distats["DiElectron"]);
  getGoodLeptonCombos(*_Muon, *_Muon, CUTS::eRMuon1, CUTS::eRMuon2, CUTS::eDiMuon, distats["DiMuon"]);

  ////Dijet cuts
  getGoodDiJets(distats["DiJet"]);

  // dit += clock() - t1;
  // t1 = clock();

  if(event % 50000 == 0) {
    cout << "Event #" << event << endl;
  }
}

////Reads cuts from Cuts.in file and see if the event has enough particles
int Analyzer::fillCuts() {
  unordered_map<string,pair<int,int> >* cut_info = histo.get_cuts();
  vector<string>* cut_order = histo.get_order();

  string cut;
  int min, max;
  bool prevTrue = true;
  int nparticles, i=0;
  int maxCut=0;

  for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
    if(isData && it->find("Gen") != string::npos) continue;
    cut = *it;
    min= cut_info->at(cut).first;
    max= cut_info->at(cut).second;
    nparticles = goodParts[ival(cut_num[cut])].size();
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      cuts_per[i]++;
      cuts_cumul[i] += (prevTrue) ? 1 : 0;
      maxCut += (prevTrue) ? 1 : 0;
    } else prevTrue = false;
  }
  return maxCut;
}


///Prints the number of events that passed each cut per event and cumulatively
//done at the end of the analysis
void Analyzer::printCuts() {
  vector<string>* cut_order = histo.get_order();
  int i =0;

  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(3);
  cout << "\n";
  cout << "Selection Efficiency " << "\n";
  cout << "Total events: " << nentries << "\n";
  cout << "         Name                     Indiv.      Cumulative\n";
  cout << "---------------------------------------------------------------------------\n";
  for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
    cout << setw(28) << *it << " ";
    if(isData && it->find("Gen") != string::npos) cout << "Skipped" << endl;
    else cout << setw(8) << cuts_per.at(i) << " (" << setw(8) << ((float)cuts_per.at(i)) / nentries << ") "
	      << setw(8) << cuts_cumul.at(i) << "( " << setw(8) << ((float)cuts_cumul.at(i)) / nentries << ") " << endl;
  }
  cout << "---------------------------------------------------------------------------\n";  

  // cout << "gen: " << gent << endl;
  // cout << "smear: " << smeart << endl;
  // cout << "trigger: " << trigt << endl;
  // cout << "reco lepton: " << recolt << endl;
  // cout << "reco1 jets: " << jet1t << endl;
  // cout << "reco2 jets: " << jet2t << endl;
  // cout << "central jets: " << cent << endl;
  // cout << "bjets: " << bjett << endl;
  // cout << "lead jet: " << leadt << endl;
  // cout << "met top: " << mett << endl;
  // cout << "susy: " << susyt << endl;
  // cout << "combo: " << combot << endl;
  // cout << "diparticle: " << dit << endl;

}

/////////////PRIVATE FUNCTIONS////////////////



///Calculates met from values from each file plus smearing and treating muons as neutrinos
void Analyzer::updateMet() {
  ////// Neutrino update before calculation
  if(distats["Run"].bmap.at("TreatMuonsAsNeutrinos")) {
    for(vec_iter it=goodParts[ival(CUTS::eRMuon1)].begin(); it!=goodParts[ival(CUTS::eRMuon1)].end(); it++) {
      if(find(goodParts[ival(CUTS::eRMuon2)].begin(), goodParts[ival(CUTS::eRMuon2)].end(), (*it)) != goodParts[ival(CUTS::eRMuon2)].end() ) continue;
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }    
    for(vec_iter it=goodParts[ival(CUTS::eRMuon2)].begin(); it!=goodParts[ival(CUTS::eRMuon2)].end(); it++) {
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }
  }
  ///---MHT and HT calculations----////
  int i=0;
  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it!=_Jet->smearP.end(); it++, i++) {
    if( (it->Pt() > distats["Run"].dmap.at("JetPtForMhtAndHt")) && (fabs(it->Eta()) < distats["Run"].dmap.at("JetEtaForMhtAndHt")) ) {
      if(distats["Run"].bmap.at("ApplyJetLooseIDforMhtAndHt") ) continue;//&& !passLooseJetID() continue;
      
      sumpxForMht -= it->Px();
      sumpyForMht -= it->Py();
      sumptForHt  += it->Pt();
    }
  }
  phiForMht = atan2(sumpyForMht,sumpxForMht);
  //  if(sumpxForMht < 0) phiForMht += (sumpyForMht >= 0) ? TMath::Pi() : -TMath::Pi();

  theMETVector.SetPxPyPzE(theMETVector.Px()+deltaMEx, theMETVector.Py()+deltaMEy, theMETVector.Pz(), 
  			  TMath::Sqrt(pow(theMETVector.Px()+deltaMEx,2) + pow(theMETVector.Py()+deltaMEy,2)));
}


/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral(TTree* BOOM, string infile) {
  BOOM->SetBranchStatus("Trigger_decision", 1);
  BOOM->SetBranchStatus("Trigger_names", 1);
  BOOM->SetBranchStatus("nTruePUInteractions", 1);
  BOOM->SetBranchStatus("bestVertices", 1);
  BOOM->SetBranchStatus("Met_type1PF_px", 1);
  BOOM->SetBranchStatus("Met_type1PF_py", 1);
  BOOM->SetBranchStatus("Met_type1PF_pz", 1);

  BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision);
  BOOM->SetBranchAddress("Trigger_names", &Trigger_names);
  BOOM->SetBranchAddress("nTruePUInteractions", &nTruePU);
  BOOM->SetBranchAddress("bestVertices", &bestVertices);
  BOOM->SetBranchAddress("Met_type1PF_px", &Met_px);
  BOOM->SetBranchAddress("Met_type1PF_py", &Met_py);
  BOOM->SetBranchAddress("Met_type1PF_pz", &Met_pz);

  read_info(FILESPACE + "ElectronTau_info.in");
  read_info(FILESPACE + "MuonTau_info.in");
  read_info(FILESPACE + "DiParticle_info.in");
  read_info(FILESPACE + "VBFCuts_info.in");
  read_info(FILESPACE + "Run_info.in");
}


///parsing method that gets info on diparts and basic run info
//put in map called "distats"
void Analyzer::read_info(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }

  vector<string> stemp;
  string group, line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      cout << "error in " << filename << "; no groups specified for data" << endl;
      exit(1);
    } else if(stemp.size() == 2) {
      if(stemp[1].find(".") != string::npos && stemp[1].find("root") == string::npos) distats[group].dmap[stemp[0]]=stod(stemp[1]);
      else if(stemp[1] == "1" || stemp[1] == "true") distats[group].bmap[stemp[0]]=true;
      else if(stemp[1] == "0" || stemp[1] == "false") distats[group].bmap[stemp[0]]=false; 
      else distats[group].smap[stemp[0]] = stemp[1];

    } else  distats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
  }
  info_file.close();
}



///Smears lepton only if specified and not a data file.  Otherwise, just filles up lorentz vectors
//of the data into the vector container smearP with is in each lepton object.
void Analyzer::smearLepton(Lepton& lepton, CUTS eGenPos, const PartStats& stats) {
  lepton.smearP.clear();

  double smearedPt;
  double smearedEta;
  double smearedPhi;
  double smearedEnergy;


  for(int i = 0; i < (int)lepton.pt->size(); i++) {
    TLorentzVector tmpSmear;
    tmpSmear.SetPtEtaPhiE(lepton.pt->at(i), lepton.eta->at(i), lepton.phi->at(i), lepton.energy->at(i));

    if(isData || !stats.bmap.at("SmearTheParticle")) {
      lepton.smearP.push_back(tmpSmear);
      continue;
    }

    TLorentzVector genVec =  matchLeptonToGen(tmpSmear, lepton.pstats["Smear"],eGenPos);
    if(genVec == TLorentzVector(0,0,0,0)) {      
      lepton.smearP.push_back(tmpSmear);
      continue;
    }

    smearedPt = (genVec.Pt()*stats.dmap.at("PtScaleOffset")) + (tmpSmear.Pt() - genVec.Pt())*stats.dmap.at("PtSigmaOffset");
    smearedEta =(genVec.Eta()*stats.dmap.at("EtaScaleOffset")) + (tmpSmear.Eta() - genVec.Eta())*stats.dmap.at("EtaSigmaOffset");
    smearedPhi = (genVec.Phi() * stats.dmap.at("PhiScaleOffset")) + (tmpSmear.Phi() - genVec.Phi())*stats.dmap.at("PhiSigmaOffset");
    smearedEnergy = (genVec.Energy()*stats.dmap.at("EnergyScaleOffset")) + (tmpSmear.Energy() - genVec.Energy())*stats.dmap.at("EnergySigmaOffset");
    
    TLorentzVector final;
    final.SetPtEtaPhiE(smearedPt, smearedEta, smearedPhi, smearedEnergy);
    lepton.smearP.push_back(final);
    deltaMEx += tmpSmear.Px() - final.Px();
    deltaMEy += tmpSmear.Py() - final.Py();

  }
}

///Same as smearlepton, just jet specific
void Analyzer::smearJet(const PartStats& stats) {
  _Jet->smearP.clear();
  TLorentzVector jetV;

  for(int i=0; i< (int)_Jet->pt->size(); i++) {
    jetV.SetPtEtaPhiE(_Jet->pt->at(i), _Jet->eta->at(i), _Jet->phi->at(i), _Jet->energy->at(i));

    if(isData || !stats.bmap.at("SmearTheJet")) {
      _Jet->smearP.push_back(jetV);
      continue;
    }
    
    if(JetMatchesLepton(*_Muon, jetV, stats.dmap.at("MuonMatchingDeltaR"), CUTS::eGMuon) ||
       JetMatchesLepton(*_Tau, jetV, stats.dmap.at("TauMatchingDeltaR"), CUTS::eGTau) ||       
       JetMatchesLepton(*_Electron, jetV,stats.dmap.at("ElectronMatchingDeltaR"), CUTS::eGElec)){

      _Jet->smearP.push_back(jetV);
      continue;
    }
    
    _Jet->smearP.push_back(stats.dmap.at("JetEnergyScaleOffset") * jetV);

  }
}

/////checks if jet is close to a lepton and the lepton is a gen particle, then the jet is a lepton object, so
//this jet isn't smeared
bool Analyzer::JetMatchesLepton(Lepton& lepton, const TLorentzVector& jetV, double partDeltaR, CUTS eGenPos) {
  TLorentzVector tempV;
  for(int j = 0; j < (int)lepton.pt->size(); j++) {
    tempV.SetPtEtaPhiE(lepton.pt->at(j), lepton.eta->at(j), lepton.phi->at(j), lepton.energy->at(j));
    if(jetV.DeltaR(tempV) < partDeltaR && matchLeptonToGen(tempV, lepton.pstats["Smear"], eGenPos) != TLorentzVector(0,0,0,0)) return true;
  }
  return false;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchLeptonToGen(const TLorentzVector& lvec, const PartStats& stats, CUTS ePos) {
  if(ePos == CUTS::eGTau) {
    return matchTauToGen(lvec, stats.dmap.at("GenMatchingDeltaR"));
  }
  TLorentzVector genVec = TLorentzVector(0,0,0,0);
  
  for(vec_iter it=goodParts[ival(ePos)].begin(); it !=goodParts[ival(ePos)].end();it++) {
    genVec.SetPtEtaPhiE(_Gen->pt->at(*it), _Gen->eta->at(*it), _Gen->phi->at(*it), _Gen->energy->at(*it));
    if(lvec.DeltaR(genVec) <= stats.dmap.at("GenMatchingDeltaR")) {
      unordered_map<string,bool>::const_iterator mother = stats.bmap.find("UseMotherID");
      if(mother != stats.bmap.end() && mother->second && abs(_Gen->motherpdg_id->at(*it)) != stats.dmap.at("MotherID")) continue; 
      return genVec;
    }
  }
  
  return TLorentzVector(0,0,0,0);
}


///Tau specific matching fucntion.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically
TLorentzVector Analyzer::matchTauToGen(const TLorentzVector& lvec, double lDeltaR) {
  TLorentzVector genVec(0,0,0,0);
  bool leptonicDecay = false;
  
  for(vec_iter it=goodParts[ival(CUTS::eGTau)].begin(); it !=goodParts[ival(CUTS::eGTau)].end();it++) {
    leptonicDecay = false;
    for(int j = 0; j < (int)_Gen->pt->size(); j++) {
      if( ((abs(_Gen->pdg_id->at(j)) == 12) || (abs(_Gen->pdg_id->at(j)) == 14)) && (_Gen->BmotherIndex->at(j) == (*it)) ) {
	leptonicDecay = true; 
	break;
      }
    }
    if(leptonicDecay) continue;

    for(vec_iter inu=goodParts[ival(CUTS::eNuTau)].begin(); inu !=goodParts[ival(CUTS::eNuTau)].end();inu++) {
      if(_Gen->BmotherIndex->at(*inu) != (*it)) continue;
      
      genVec.SetPtEtaPhiE(_Gen->pt->at(*it)-_Gen->pt->at(*inu), _Gen->eta->at(*it)-_Gen->eta->at(*inu), 
			   _Gen->phi->at(*it)-_Gen->phi->at(*inu), _Gen->energy->at(*it)-_Gen->energy->at(*inu));
      if(lvec.DeltaR(genVec) <= lDeltaR) {
	return genVec;
      }
    }
  }

  return TLorentzVector(0,0,0,0);
}


////Calculates the number of gen particles.  Based on id number and status of each particle
void Analyzer::getGoodGen(int particle_id, int particle_status, CUTS ePos, const PartStats& stats) {
  for(int j = 0; j < (int)_Gen->pt->size(); j++) {
    if(particle_id == 15 && (_Gen->pt->at(j) < stats.dmap.at("TauPtMinCut")) && (abs(_Gen->eta->at(j)) > stats.dmap.at("TauEtaMaxCut"))) continue;
    
    if((abs(_Gen->pdg_id->at(j)) == particle_id) && (_Gen->status->at(j) == particle_status)) {
      goodParts[ival(ePos)].push_back(j);
    }
  }
}

////Tau neutrino specific function used for calculating the number of hadronic taus
void Analyzer::getGoodTauNu() {
  for(int j = 0; j < (int)_Gen->pt->size(); j++) {
    
    if( (abs(_Gen->pdg_id->at(j)) == 16) && (abs(_Gen->pdg_id->at(_Gen->BmotherIndex->at(j))) == 15) && (_Gen->status->at(_Gen->BmotherIndex->at(j)) == 2) ) {
      goodParts[ival(CUTS::eNuTau)].push_back(j);

    }
  }
}

///Function used to find the number of reco leptons that pass the various cuts.
///Divided into if blocks for the different lepton requirements.
void Analyzer::getGoodRecoLeptons(Lepton& lep, CUTS ePos, CUTS eGenPos, const PartStats& stats) {
  int i = 0;

  for(vector<TLorentzVector>::iterator it=lep.smearP.begin(); it != lep.smearP.end(); it++, i++) {
    TLorentzVector lvec = (*it);
    
    if (fabs(lvec.Eta()) > stats.dmap.at("EtaCut")) continue;
    if (lvec.Pt() < stats.pmap.at("PtCut").first || lvec.Pt() > stats.pmap.at("PtCut").second) continue;

    if((lep.pstats.at("Smear").bmap.at("MatchToGen")) && (!isData)) {   /////check
      if(matchLeptonToGen(lvec, lep.pstats["Smear"] ,eGenPos) == TLorentzVector(0,0,0,0)) continue;
    }

    if(ePos == CUTS::eRMuon1 || ePos == CUTS::eRMuon2) {      ////////////////MUON CUTS/////////////
      Muon& partM = static_cast<Muon&>(lep);
    
      if(stats.bmap.at("DoDiscrByTightID") && (partM.tight->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrBySoftID") && (partM.soft->at(i) == 0)) continue;
      
      if (stats.bmap.at("DoDiscrByIsolation")) {
	double maxIsoval = max(0.0, partM.isoNeutralHadron->at(i) + partM.isoPhoton->at(i) - 0.5 * partM.isoPU->at(i));
	double isoSum = (partM.isoCharged->at(i) + maxIsoval) / lvec.Pt();
	if(isoSum < stats.pmap.at("IsoSumPtCutValue").first || isoSum >= stats.pmap.at("IsoSumPtCutValue").second) continue;
      }
    } else if(ePos == CUTS::eRElec1 || ePos == CUTS::eRElec2) {    ///////////////ELECTRON CUT///////////
      Electron& partE = static_cast<Electron&>(lep);

      //----Require electron to pass ID discriminators
      if(stats.bmap.at("DoDiscrByVetoID") && (partE.isPassVeto->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByLooseID") && (partE.isPassLoose->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByMediumID") && (partE.isPassMedium->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByTightID") && (partE.isPassTight->at(i) == 0)) continue;
      if(stats.bmap.at("DoDiscrByHEEPID") && (partE.isPassHEEPId->at(i) == 0)) continue;

      if (stats.bmap.at("DoDiscrByIsolation")) {
	double maxIsoval = std::max(0.0 , partE.isoNeutralHadrons->at(i) + partE.isoPhotons->at(i) - 0.5 * partE.isoPU->at(i) );
	double isoSum = (partE.isoChargedHadrons->at(i) + maxIsoval) / lvec.Pt();
	if(isoSum < stats.pmap.at("IsoSumPtCutValue").first || isoSum >= stats.pmap.at("IsoSumPtCutValue").second) continue;
      }

    } else if(ePos == CUTS::eRTau1 || ePos == CUTS::eRTau2) {   /////////////TAU CUT/////////////////
      Taus& partT = static_cast<Taus&>(lep);
      if (stats.bmap.at("DoDiscrByLeadTrack")) {
	if(partT.leadChargedCandPt->at(i) < stats.dmap.at("LeadTrackThreshold")) continue;
      }

      // ----Isolation requirement
      if (stats.bmap.at("DoDiscrByIsolation")) {
	//--- max isolation requirement
	maxIso = (ePos == CUTS::eRTau1) ? partT.maxIso.first->at(i) : partT.maxIso.second->at(i);
	if(maxIso < 0.5) continue;
	if(stats.smap.at("DiscrByMinIsolation") != "ZERO") {
	  minIso = (ePos == CUTS::eRTau1) ? partT.minIso.first->at(i) : partT.minIso.second->at(i);
	  if (minIso > 0.5) continue;
	}
      }

      // ----Require 1 or 3 prongs
      if(stats.smap.at("DiscrByProngType").find("hps") != string::npos && partT.decayModeFindingNewDMs->at(i) < 0.5) continue;
      if(!passProng(stats.smap.at("DiscrByProngType"), partT.nProngs->at(i))) continue;

      // ----Electron and Muon vetos
      if (stats.bmap.at("DoDiscrAgainstElectron") && !stats.bmap.at("SelectTausThatAreElectrons")) {
	againstElectron = (ePos == CUTS::eRTau1) ? partT.againstElectron.first->at(i) : partT.againstElectron.second->at(i);
	if (againstElectron < 0.5) continue;
      } else if (stats.bmap.at("SelectTausThatAreElectrons")) {
	if (againstElectron > 0.5) continue;
      }

      if (stats.bmap.at("DoDiscrAgainstMuon") && !stats.bmap.at("SelectTausThatAreMuons")) {
	againstMuon = (ePos == CUTS::eRTau1) ? partT.againstMuon.first->at(i) : partT.againstMuon.second->at(i);
	if (againstMuon < 0.5) continue;
      } else if (stats.bmap.at("SelectTausThatAreMuons")) {
	againstMuon = (ePos == CUTS::eRTau1) ? partT.againstMuon.first->at(i) : partT.againstMuon.second->at(i);
	if (againstMuon > 0.5) continue;
      }

      if (stats.bmap.at("DoDiscrByCrackCut") && isInTheCracks(lvec.Eta())) continue;

      // ----anti-overlap requirements
      if (stats.bmap.at("RemoveOverlapWithMuon1s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithMuon2s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithElectron1s") && isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"))) continue;
      if (stats.bmap.at("RemoveOverlapWithElectron2s") && isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"))) continue;
    }
    goodParts[ival(ePos)].push_back(i);    
  }

}

////Jet specific function for finding the number of jets that pass the cuts.
//used to find the nubmer of good jet1, jet2, central jet, 1st and 2nd leading jets and bjet.
void Analyzer::getGoodRecoJets(CUTS ePos, const PartStats& stats) {
  int i=0;

  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it != _Jet->smearP.end(); it++, i++) {
    TLorentzVector lvec = (*it);
    ///if else loop for central jet requirements

    if( ePos == CUTS::eRCenJet) {
      if(lvec.Pt() > 2.5) continue;
    } else if (fabs(lvec.Eta()) < stats.pmap.at("EtaCut").first || fabs(lvec.Eta()) > stats.pmap.at("EtaCut").second) continue;
    if (lvec.Pt() < stats.dmap.at("PtCut")) continue;

    /// BJet specific
    if(ePos == CUTS::eRBJet) {
      if(stats.bmap.at("ApplyJetBTagging") && _Jet->bDiscriminator->at(i) <= stats.dmap.at("JetBTaggingCut")) continue;
      if((stats.bmap.at("MatchBToGen")) && !isData && abs(_Jet->partonFlavour->at(i)) != 5) continue;
    } else if (stats.bmap.at("ApplyLooseID") && !passedLooseJetID(i)) continue;

  // ----anti-overlap requirements
    if(stats.bmap.at("RemoveOverlapWithMuon1s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithMuon2s") && isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithElectron1s") && isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithElectron2s") && isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithTau1s") && isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"))) continue;
    if(stats.bmap.at("RemoveOverlapWithTau2s") && isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"))) continue;

    /////fill up array
    goodParts[ival(ePos)].push_back(i);    
  }
  
  if(ePos == CUTS::eR1stJet || ePos == CUTS::eR2ndJet) {
    int potential = -1;
    double prevPt = -1;
    for(vec_iter leadit = goodParts[ival(ePos)].begin(); leadit != goodParts[ival(ePos)].end(); ++leadit) {
      if(((ePos == CUTS::eR2ndJet && (*leadit) != leadIndex) || ePos == CUTS::eR1stJet) && _Jet->smearP.at(*leadit).Pt() > prevPt) {
	potential = (*leadit);
	prevPt = _Jet->smearP.at(*leadit).Pt();
      }
    }
    goodParts[ival(ePos)].clear();
    goodParts[ival(ePos)].push_back(potential);
  }

} 

///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(vec_iter it=goodParts[ival(ePos)].begin(); it < goodParts[ival(ePos)].end(); it++) {
    if(lvec.DeltaR(overlapper.smearP.at(*it)) < MatchingDeltaR) return true;
  }
  return false;
}

///Tests if tau decays into the specified number of jet prongs.
bool Analyzer::passProng(string prong, int value) {
  return ( (prong.find("1") != string::npos && value == 1) ||
	   (prong.find("2") != string::npos && value == 2) ||
	   (prong.find("3") != string::npos && value == 3) );
}


////Tests if tau is within the cracks of the detector (the specified eta ranges)
bool Analyzer::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 ||
          (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
          (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
          (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
          (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}
 

//Tests if a jet meets a litany of different tests
bool Analyzer::passedLooseJetID(int nobj) {
  if (_Jet->neutralHadEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->neutralEmEmEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->numberOfConstituents->at(nobj) <= 1) return false;
  if (_Jet->muonEnergyFraction->at(nobj) >= 0.80) return false;
  if ( (fabs(_Jet->smearP.at(nobj).Eta()) < 2.4) && 
           ((_Jet->chargedHadronEnergyFraction->at(nobj) <= 0.0) || 
	    (_Jet->chargedMultiplicity->at(nobj) <= 0.0) || 
	    (_Jet->chargedEmEnergyFraction->at(nobj) >= 0.99) )) return false;
  return true;
}


///sees if the event passed one of the two cuts provided
bool Analyzer::passTriggerCuts(string TriggerN) {
  if(prevTrig[TriggerN].first >= (int)Trigger_names->size() ||
     distats["Run"].smap[TriggerN+"FirstRequirement"] != Trigger_names->at(prevTrig[TriggerN].first)) {

    int temp = 0;
    vector<string>::iterator it = find(Trigger_names->begin(), Trigger_names->end(), distats["Run"].smap[TriggerN+"FirstRequirement"]);
    while(it != Trigger_names->begin()) {
      temp++; it--;
    }
    prevTrig[TriggerN].first = temp;
  }
  if(prevTrig[TriggerN].second >= (int)Trigger_names->size() ||
     distats["Run"].smap[TriggerN+"SecondRequirement"] != Trigger_names->at(prevTrig[TriggerN].second)) {
    int temp = 0;
    vector<string>::iterator it = find(Trigger_names->begin(), Trigger_names->end(), distats["Run"].smap[TriggerN+"SecondRequirement"]);
    while(it != Trigger_names->begin()) {
      temp++; 
      it--;
    }
    prevTrig[TriggerN].second = temp;
  }     
  if( (prevTrig[TriggerN].first < (int)Trigger_names->size() && Trigger_decision->at(prevTrig[TriggerN].first) == 1) || 
	(prevTrig[TriggerN].second < (int)Trigger_names->size() && Trigger_decision->at(prevTrig[TriggerN].second) == 1) ) return true;

  return false;
}

////VBF specific cuts dealing with the leading jets.
void Analyzer::VBFTopologyCut() {
  PartStats stats = distats["VBFSUSY"];
  TLorentzVector ljet1 = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0));
  TLorentzVector ljet2 = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));
  
  if(stats.bmap.at("DiscrByMass")) {
    if((ljet1 + ljet2).M() < stats.pmap.at("MassCut").first) return;
    if((ljet1 + ljet2).M() > stats.pmap.at("MassCut").second) return;
  }

  if(stats.bmap.at("DiscrByPt")) {
    if((ljet1 + ljet2).Pt() < stats.pmap.at("PtCut").first) return;
    if((ljet1 + ljet2).Pt() > stats.pmap.at("PtCut").second) return;
  }
  
  if(stats.bmap.at("DiscrByDeltaEta")) {
    if(fabs(ljet1.Eta() - ljet2.Eta()) < stats.pmap.at("DeltaEtaCut").first) return;
    if(fabs(ljet1.Eta() - ljet2.Eta()) > stats.pmap.at("DeltaEtaCut").second) return;
  }

  if(stats.bmap.at("DiscrByDeltaPhi")) {
    if(absnormPhi(ljet1.Phi() - ljet2.Phi()) < stats.pmap.at("DeltaPhiCut").first) return;
    if(absnormPhi(ljet1.Phi() - ljet2.Phi()) > stats.pmap.at("DeltaPhiCut").second) return;
  }
  
  if(stats.bmap.at("DiscrByOSEta")) {
    if((ljet1.Eta() * ljet2.Eta()) >= 0) return;
  }
  if(stats.bmap.at("DiscrByMet")) {
    if(theMETVector.Pt() < stats.pmap.at("RecoMetCut").first) return;
    if(theMETVector.Pt() > stats.pmap.at("RecoMetCut").second) return;
  }
  if(stats.bmap.at("DiscrByMHT") && sqrt(pow(sumpxForMht,2) + pow(sumpyForMht,2)) < stats.dmap.at("MhtCut")) return;

  if(stats.bmap.at("DiscrByHT") && sumptForHt < stats.dmap.at("HtCut")) return; 

  double dphi1 = normPhi(ljet1.Phi() - theMETVector.Phi());
  double dphi2 = normPhi(ljet2.Phi() - theMETVector.Phi());
  double r1, r2, alpha;
  
  if(stats.bmap.at("DiscrByR1")) {
    r1 = sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) );
    if(r1 < stats.pmap.at("R1Cut").first || r1 > stats.pmap.at("R1Cut").second) return;

  }
  if(stats.bmap.at("DiscrByR2")) {
    r2 = sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0) );
    if(r2 < stats.pmap.at("R2Cut").first || r2 > stats.pmap.at("R2Cut").second) return;
  }
  if(stats.bmap.at("DiscrByAlpha")) {
    TLorentzVector addVec = ljet1 + ljet2;
    alpha = (addVec.M() > 0) ? ljet2.Pt() / addVec.M() : -1;
    if(alpha < stats.pmap.at("AlphaCut").first || alpha > stats.pmap.at("AlphaCut").second) return;
  }
  if(stats.bmap.at("DiscrByDphi1")) {
    if(abs(dphi1) < stats.pmap.at("Dphi1Cut").first || abs(dphi1) > stats.pmap.at("Dphi1Cut").second) return;
  }
  if(stats.bmap.at("DiscrByDphi2")) {
    if(abs(dphi2) < stats.pmap.at("Dphi2Cut").first || abs(dphi2) > stats.pmap.at("Dphi2Cut").second) return;
  }

  goodParts[ival(CUTS::eSusyCom)].push_back(0);
}


void Analyzer::getGoodMetTopologyLepton(Lepton& lep, CUTS eReco, CUTS ePos, const PartStats& stats) {

  for(vec_iter it=goodParts[ival(eReco)].begin(); it != goodParts[ival(eReco)].end(); it++) {
    if ((ePos != CUTS::eTTau1 && ePos != CUTS::eTTau2) && stats.bmap.at("DiscrIfIsZdecay")) { 
      if(isZdecay(lep.smearP.at(*it), lep)) continue;
    }

    if (stats.bmap.at("DiscrByMetDphi")) {
      double Dphi = absnormPhi(lep.smearP.at(*it).Phi() - theMETVector.Phi());
      if(Dphi < stats.pmap.at("MetDphiCut").first || Dphi > stats.pmap.at("MetDphiCut").second) continue;
    }

    if (stats.bmap.at("DiscrByMetMt")) {
      double LeptonMetMt = calculateLeptonMetMt(lep.smearP.at(*it));
      if( LeptonMetMt < stats.pmap.at("MetMtCut").first || LeptonMetMt > stats.pmap.at("MetMtCut").second) continue;
    }
    goodParts[ival(ePos)].push_back(*it);
  }
}


//-----Calculate lepton+met transverse mass
double Analyzer::calculateLeptonMetMt(const TLorentzVector& Tobj) {
  double px = Tobj.Px() + theMETVector.Px();
  double py = Tobj.Py() + theMETVector.Py();
  double et = Tobj.Et() + theMETVector.Energy(); //TMath::Sqrt((theMETVector.Px() * theMETVector.Px()) + (theMETVector.Py() * theMETVector.Py()));
  double mt2 = et*et - (px*px + py*py);
  return (mt2 >= 0) ? sqrt(mt2) : -1;
}

/////all it does is add lorentz vectors. 
/////keep in case needed later
// TLorentzVector Analyzer::CalculateTheDiJet4Momentum(TLorentzVector* Tobj1, TLorentzVector* Tobj2) {
//   return (*Tobj1) + (*Tobj2);
// }


/////Calculate the diparticle mass based on how to calculate it
///can use Collinear Approximation, which can fail (number failed available in a histogram)
///can use VectorSumOfVisProductAndMet which is sum of particles and met
///Other which is adding without met
double Analyzer::diParticleMass(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string howCalc) {
  bool ratioNotInRange = false;
  TLorentzVector The_LorentzVect;


  //////check this equation/////
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + theMETVector.Px())) - (Tobj2.Px() * (Tobj1.Py() + theMETVector.Py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + theMETVector.Py())) - (Tobj1.Py() * (Tobj2.Px() + theMETVector.Px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    ratioNotInRange=!((x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px() / x1) + (Tobj2.Px() / x2), (Tobj1.Py() / x1) + (Tobj2.Py() / x2), (Tobj1.Pz() / x1) + (Tobj2.Pz() / x2), (Tobj1.Energy() / x1) + (Tobj2.Energy() / x2) );
      return The_LorentzVect.M();
    }
  } 

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    double px = Tobj1.Px() + Tobj2.Px() + theMETVector.Px();
    double py = Tobj1.Py() + Tobj2.Py() + theMETVector.Py();
    double pz = Tobj1.Pz() + Tobj2.Pz();
    double e  = Tobj1.Energy() + Tobj2.Energy() + theMETVector.Energy();///TMath::Sqrt(pow(theMETVector.Px(),2) + pow(theMETVector.Py(),2));
    The_LorentzVect.SetPxPyPzE(px, py, pz, e);
    return The_LorentzVect.M();
  }

  return (Tobj1 + Tobj2).M();
}

////Tests if the CollinearApproximation works for finding the mass of teh particles
bool Analyzer::passDiParticleApprox(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, string howCalc) {
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + theMETVector.Px())) - (Tobj2.Px() * (Tobj1.Py() + theMETVector.Py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + theMETVector.Py())) - (Tobj1.Py() * (Tobj2.Px() + theMETVector.Px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    return (x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.);
  } else {
    return true;
  }
}




/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonCombos(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats) {
  for(vec_iter i1=goodParts[ival(ePos1)].begin(); i1 != goodParts[ival(ePos1)].end(); i1++) {
    for(vec_iter i2=goodParts[ival(ePos2)].begin(); i2 != goodParts[ival(ePos2)].end(); i2++) {
      if(stats.bmap.at("DiscrByDeltaR") && (lep1.smearP.at(*i1).DeltaR(lep2.smearP.at(*i2))) < stats.dmap.at("DeltaRCut")) continue;
   
      if(stats.smap.at("DiscrByOSLSType") == "LS" && (lep1.charge->at(*i1) * lep2.charge->at(*i2) >= 0)) continue;
      else if(stats.smap.at("DiscrByOSLSType") == "OS" && (lep1.charge->at(*i1) * lep2.charge->at(*i2) <= 0)) continue;

      if(stats.bmap.at("DiscrByCosDphi")) {
	double Dphi = absnormPhi( lep1.smearP.at(*i1).Phi() - lep2.smearP.at(*i2).Phi());
	if(cos(Dphi) > stats.pmap.at("CosDphiCut").first || cos(Dphi) < stats.pmap.at("CosDphiCut").second) continue;
      }
  // ----Mass window requirement

      if (stats.bmap.at("DiscrByMassReco")) {
      	double diMass = diParticleMass(lep1.smearP.at(*i1),lep2.smearP.at(*i2), stats.smap.at("HowCalculateMassReco"));
      	if( diMass < stats.pmap.at("MassCut").first || diMass > stats.pmap.at("MassMaxCut").second) continue;
      }

      if (stats.bmap.at("DiscrByCDFzeta2D")) {
      	double CDFzeta = stats.dmap.at("PZetaCutCoefficient") * getPZeta(lep1.smearP.at(*i1), lep2.smearP.at(*i2)) 
	  + stats.dmap.at("PZetaVisCutCoeffiecient") * getPZetaVis(lep1.smearP.at(*i1), lep2.smearP.at(*i2));
      	if( CDFzeta < stats.pmap.at("CDFzeta2DCutValue").first || CDFzeta > stats.pmap.at("CDFzeta2DCutValue").second ) continue;
      }
      //////////abs on the difference????
      ///////////////////
      if (stats.bmap.at("DiscrByDeltaPtDivSumPt")) {
	double ptDiv = (lep1.smearP.at(*i1).Pt() - lep2.smearP.at(*i2).Pt()) / (lep1.smearP.at(*i1).Pt() + lep2.smearP.at(*i2).Pt());
	if( ptDiv < stats.pmap.at("DeltaPtDivSumPtCut").first || ptDiv > stats.pmap.at("DeltaPtDivSumPtCut").second) continue;
      }

      if (stats.bmap.at("DiscrByDeltaPt")) {
	double deltaPt = lep1.smearP.at(*i1).Pt() - lep2.smearP.at(*i2).Pt(); 
	if(deltaPt < stats.pmap.at("DeltaPtCutValue").first || deltaPt > stats.pmap.at("DeltaPtCutValue").second) continue;
      }
      ///Particlesp that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2 
      goodParts[ival(ePosFin)].push_back((*i1)*_Gen->pt->size() + (*i2));
    }
  }
}


/////Same as gooddilepton, just jet specific
void Analyzer::getGoodDiJets(const PartStats& stats) {
  // ----Separation cut between jets (remove overlaps)
  for(vec_iter ij1=goodParts[ival(CUTS::eRJet1)].begin(); ij1 != goodParts[ival(CUTS::eRJet1)].end(); ij1++) {
    for(vec_iter ij2=goodParts[ival(CUTS::eRJet2)].begin(); ij2 != goodParts[ival(CUTS::eRJet2)].end(); ij2++) {

      if (stats.bmap.at("DiscrByDeltaR")) {
	if(_Jet->smearP.at(*ij1).DeltaR(_Jet->smearP.at(*ij2)) < stats.dmap.at("DeltaRCut")) continue;
      }

      if (stats.bmap.at("DiscrByDeltaEta")) {
	if(fabs(_Jet->smearP.at(*ij1).Eta() - _Jet->smearP.at(*ij2).Eta()) < stats.pmap.at("DeltaEtaCut").first) continue;
	if(fabs(_Jet->smearP.at(*ij1).Eta() - _Jet->smearP.at(*ij2).Eta()) > stats.pmap.at("DeltaEtaCut").second) continue;
      }

      if (stats.bmap.at("DiscrByDeltaPhi")) {
	if(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi()) < stats.pmap.at("DeltaPhiCut").first) continue;
	if(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi()) > stats.pmap.at("DeltaPhiCut").second) continue;
      } 

      if (stats.bmap.at("DiscrByOSEta")) {
	if((_Jet->smearP.at(*ij1).Eta() * _Jet->smearP.at(*ij2).Eta()) >= 0) continue;
      }
      // ----Require both legs to be almost back-to-back in phi
      if (stats.bmap.at("DiscrByCosDphi")) {
	if(cos(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi())) > stats.pmap.at("CosDphiMaxCut").first) continue;
	if(cos(absnormPhi(_Jet->smearP.at(*ij1).Phi() - _Jet->smearP.at(*ij2).Phi())) < stats.pmap.at("CosDphiMinCut").second) continue;
      }
      // ----Mass window requirement
      if (stats.bmap.at("DiscrByMassReco")) {
	if( ((_Jet->smearP.at(*ij1) + _Jet->smearP.at(*ij2)).M() < stats.pmap.at("MassCut").first) || ((_Jet->smearP.at(*ij1) + _Jet->smearP.at(*ij2)).M() > stats.pmap.at("MassCut").second) ) continue;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2 
      goodParts[ival(CUTS::eDiJet)].push_back((*ij1)*_Jet->smearP.size() + (*ij2));
    }
  }
}


///////Only tested for if is Zdecay, can include massptasymmpair later?
/////Tests to see if a light lepton came form a zdecay
bool Analyzer::isZdecay(const TLorentzVector& theObject, Lepton& lep) {
  bool eventIsZdecay = false;
  const float zMass = 90.1876;
  const float zWidth = 2.4952;
  float zmmPtAsymmetry = -10.;

  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.
  for(vector<TLorentzVector>::iterator lepit= lep.smearP.begin(); lepit != lep.smearP.end(); lepit++) {
    if(theObject.DeltaR(*lepit) < 0.3) continue;
    if(theObject == (*lepit)) continue;

    TLorentzVector The_LorentzVect = theObject + (*lepit);
    zmmPtAsymmetry = (theObject.Pt() - lepit->Pt()) / (theObject.Pt() + lepit->Pt());
    //    theMassPtAsymmPair = std::make_pair<float, float>(The_LorentzVect.M(), float(zmmPtAsymmetry));

    if( (abs(The_LorentzVect.M() - zMass) < 3.*zWidth) || (fabs(zmmPtAsymmetry) < 0.20) ) {
      eventIsZdecay = true;
      break;
    }
  }

  return eventIsZdecay;
}


///Calculates the Pzeta value
double Analyzer::getPZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double px = visPx + theMETVector.Px();
  double py = visPy + theMETVector.Py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}


///calculates the visible pzeta value
double Analyzer::getPZetaVis(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}

///Normalizes phi to be between -PI and PI
double Analyzer::normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}

///Takes the absolute value of of normPhi (made because constant use)
double Analyzer::absnormPhi(double phi) {
  return abs(normPhi(phi));
}



// void BSM3GAnalyzer::fillHistograms(string group, vector<int> group_list) {
//   for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++){
//     double weight = isrgluon_weight * pdfWeightVector.at(NpdfID);
//       _hNPVertices[i][NpdfID]->Fill(bestVertices,weight);
//     }

    
    
//       int nGenTaus = 0;
//       for(int j = 0; j < Gen_pt->size(); j++) {
//         if(abs(Gen_pdg_id->at(j)) == 15) {

//           _hGenTauStatusCode[i][NpdfID]->Fill(Gen_status->at(j),weight);
//         }
//         if((abs(Gen_pdg_id->at(j)) == 15) && (Gen_status->at(j) == 2)) {
//           nGenTaus++;
// 	}
//       }
    
//       _hNGenTau[i][NpdfID]->Fill(nGenTaus,weight);

//       int nGenHadTaus = 0;
//       bool IsItAHadronicDecay; 
//       TLorentzVector theGenObject(0,0,0,0); 
//       TLorentzVector theNeutrinoObject(0,0,0,0);
//       vector<bool> IsItAHadronicDecayVector; 
//       IsItAHadronicDecayVector.clear();
//       vector<int> tempTauIndexVector; 
//       tempTauIndexVector.clear();
//       vector<TLorentzVector> tempNeutrinoMomentumVector; 
//       tempNeutrinoMomentumVector.clear();
//       for(int j = 0; j < Gen_pt->size(); j++) {
//         if( (abs(Gen_pdg_id->at(j)) == 16) && (abs(Gen_pdg_id->at(Gen_BmotherIndex->at(j))) == 15) && (Gen_status->at(Gen_BmotherIndex->at(j)) == 2) ) {
//           tempTauIndexVector.push_back(Gen_BmotherIndex->at(j));
//           theNeutrinoObject.SetPtEtaPhiE(Gen_pt->at(j), Gen_eta->at(j), Gen_phi->at(j), Gen_energy->at(j));
//           tempNeutrinoMomentumVector.push_back(theNeutrinoObject);
//         }
//       }
//       if(tempTauIndexVector.size() > 0) {
//         for(int jj = 0; jj < tempTauIndexVector.size(); jj++) {
//           IsItAHadronicDecay = true;
//           for(int j = 0; j < Gen_pt->size(); j++) {
//             if( ((abs(Gen_pdg_id->at(j)) == 12) || (abs(Gen_pdg_id->at(j)) == 14)) && (Gen_BmotherIndex->at(j) == tempTauIndexVector.at(jj)) ) {
//               IsItAHadronicDecay = false; // it is not a hadronic tau decay since it decayed to a electron/muon neutrino
//             }
//           }
//           IsItAHadronicDecayVector.push_back(IsItAHadronicDecay);
//         }
//         for(int jj = 0; jj < tempTauIndexVector.size(); jj++) {
//           for(int j = 0; j < Gen_pt->size(); j++) {
//             if(j == tempTauIndexVector.at(jj)) {
//               theGenObject.SetPtEtaPhiE(Gen_pt->at(j), Gen_eta->at(j), Gen_phi->at(j), Gen_energy->at(j));
//               theGenObject = theGenObject - tempNeutrinoMomentumVector.at(jj);
//               if( (IsItAHadronicDecayVector.at(jj)) ) {
//                 nGenHadTaus++;
//                 _hGenHadTauPt[i][NpdfID]->Fill(theGenObject.Pt(),weight);
//                 _hGenHadTauEta[i][NpdfID]->Fill(theGenObject.Eta(),weight);
//               }
//             }
//           }
//         }
//       }
//       _hNGenHadTau[i][NpdfID]->Fill(nGenHadTaus,weight);

//       _hNGenMuon[i][NpdfID]->Fill(nGenMuons,weight);

//       for(int j = 0; j < Gen_pt->size(); j++) {
//         if(abs(Gen_pdg_id->at(j)) == 32) {
//             _hGenZprimeStatusCode[i][NpdfID]->Fill(Gen_status->at(j),weight);
//         }
//         if((abs(Gen_pdg_id->at(j)) == 32) && (Gen_status->at(j) != 22)) {
//           TLorentzVector genObjt1;
//           genObjt1.SetPtEtaPhiE(Gen_pt->at(j), Gen_eta->at(j), Gen_phi->at(j), Gen_energy->at(j));
//           _hGenZprimeMass[i][NpdfID]->Fill(genObjt1.M(),weight);
//         }
//       }


//       for(int j = 0; j < Gen_pt->size(); j++) {
//         for(int jj = 0; jj < Gen_pt->size(); jj++) {
//           if((abs(Gen_pdg_id->at(j)) == 15) && (Gen_status->at(j) == 2) && (abs(Gen_pdg_id->at(jj)) == 15) && (Gen_status->at(jj) == 2) && (jj > j)) {
//             TLorentzVector genObjt1;
//             genObjt1.SetPtEtaPhiE(Gen_pt->at(j), Gen_eta->at(j), Gen_phi->at(j), Gen_energy->at(j));
//             TLorentzVector genObjt2;
//             genObjt2.SetPtEtaPhiE(Gen_pt->at(jj), Gen_eta->at(jj), Gen_phi->at(jj), Gen_energy->at(jj));
//             _hGenDiTauMass[i][NpdfID]->Fill(CalculateTheDiJet4Momentum(genObjt1,genObjt2).second.M(),weight);
//           }
//         }
//       }

//     }



//       // ------Central Jet Histograms
//       int nCJets = 0;
//       for(int j = 0; j < Jet_pt->size(); j++) {
// 	if (!passRecoCentralJetCuts(j)) continue;
// 	_hCentralJetPt[i][NpdfID]->Fill(smearedJetMomentumVector.at(j).Pt(),weight);
// 	_hCentralJetEta[i][NpdfID]->Fill(smearedJetMomentumVector.at(j).Eta(),weight);
// 	nCJets++;

//       _hNCentralJet[i][NpdfID]->Fill(nCJets,weight);

  
  // if(group == "FillGenTau") {
  //   for(int i = 0; i < goodParts[ival(eGTau)].size(); i++) {
  //     histo.addVal(Gen->energy->at(i), group, );   //GenTauEnergy
  //     histo.addVal(Gen->pt->at(i), group, );   //GenTauPt
  //     histo.addVal(Gen->eta->at(i),group, );   //GenTauEta
  //     histo.addVal(Gen->phi->at(i),group, );   //GenTauPhi
  //   }

  // } else if(group == "FillGenMuon") {
  //   for(int i = 0; i < goodParts[ival(eGMuon)].size(); i++) {

  // 	switch(j) {

  // 	case 0: Fill(Gen->energy->at(i));   //GenMuonEnergy
  //       case 1: Fill(Gen_pt->at(i));   //GenMuonPt
  //       case 2: Fill(Gen_eta->at(i));   //GenMuonEta
  //       case 3: Fill(Gen_phi->at(i));   //GenMuonPhi

  // 	}
  //     }
  //   }

  //     if(Muon->smearP.at(j).Pt() >= leadingmuonpt) {
  // 	leadingmuonpt = Muon->smearP.at(j).Pt();
  // 	leadingmuoneta = Muon->smearP.at(j).Eta();
  //   if(nMuons > 0) {
  //     _hFirstLeadingMuon1Pt->Fill(leadingmuonpt);
  //     _hFirstLeadingMuon1Eta->Fill(leadingmuoneta);
  //   }
  //////////num with pt/eta/phi/e 

////Grabs a list of the groups of histograms to be filled and asked Fill_folder to fill up the histograms
void Analyzer::fill_histogram() {
  int maxCut = fillCuts();
  vector<string> groups = *histo.get_groups();
  wgt = pu_weight;
  for(vector<string>::iterator it = groups.begin(); it!=groups.end(); it++) {
    fill_Folder(*it, maxCut);
  }
}

///Function that fills up the histograms
void Analyzer::fill_Folder(string group, int max) {

  if(group == "FillTauJet1" || group == "FillTauJet2" || group == "FillMuon1" || group == "FillMuon2" || group == "FillJet1" || group == "FillJet2" || group == "FillBJet") {
    Particle* part;
    if(group == "FillTauJet1" || group == "FillTauJet2") part=_Tau;
    else if(group == "FillMuon1" || group == "FillMuon2") part=_Muon;
    else part = _Jet;
    CUTS ePos = fill_num[group];
    //cout << goodParts[ival(ePos)].size() << " " << goodParts[ival(CUTS::eRTau1)].size() << endl;
    for(vec_iter it=goodParts[ival(ePos)].begin(); it!=goodParts[ival(ePos)].end(); it++) {
      histo.addVal(part->smearP.at(*it).Energy(), group,max, "Energy", wgt);
      histo.addVal(part->smearP.at(*it).Pt(), group,max, "Pt", wgt);
      histo.addVal(part->smearP.at(*it).Eta(), group,max, "Eta", wgt);
      histo.addVal(part->smearP.at(*it).Phi(), group,max, "Phi", wgt);
      if(dynamic_cast<Taus*>(part) != NULL) {
	histo.addVal(_Tau->nProngs->at(*it), group,max, "NumSignalTracks", wgt);
  	histo.addVal(_Tau->charge->at(*it), group,max, "Charge", wgt);
	histo.addVal(_Tau->leadChargedCandPt->at(*it), group,max, "SeedTracks", wgt);
      } else if(dynamic_cast<Muon*>(part)) {
	histo.addVal(calculateLeptonMetMt(_Muon->smearP.at(*it)), group,max, "MetMt", wgt);  
      }
    }

    histo.addVal(goodParts[ival(ePos)].size(), group,max, "N", wgt);

  } else if(group == "FillSusyCuts") {

    histo.addVal(sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)),group,max, "MHT", wgt);
    histo.addVal(sumptForHt,group,max, "HT", wgt);  
    histo.addVal(sumptForHt + sqrt((sumpxForMht * sumpxForMht) + (sumpyForMht * sumpyForMht)),group,max, "Meff", wgt);
    histo.addVal(theMETVector.Pt(), group,max, "Met", wgt);
    if(goodParts[ival(CUTS::eR1stJet)].at(0) !=-1 && goodParts[ival(CUTS::eR2ndJet)].at(0) != -1) {
      TLorentzVector DiJet = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0)) + _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));
      histo.addVal(absnormPhi(theMETVector.Phi() - DiJet.Phi()), group,max, "MetDiJetDeltaPhi", wgt);
    }
  

  } else if(group == "FillLeadingJet" && goodParts[ival(CUTS::eSusyCom)].size() != 0) {
    TLorentzVector first = _Jet->smearP.at(goodParts[ival(CUTS::eR1stJet)].at(0));
    TLorentzVector second = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0));

    histo.addVal(first.Pt(),group,max, "FirstPt", wgt);
    histo.addVal(second.Pt(),group,max, "SecondPt", wgt);
    histo.addVal(first.Eta(),group,max, "FirstEta", wgt);
    histo.addVal(second.Eta(),group,max, "SecondEta", wgt);
    
    TLorentzVector LeadDiJet = first + second;
    
    histo.addVal(LeadDiJet.M(), group,max, "Mass", wgt); 
    histo.addVal(LeadDiJet.Pt(), group,max, "Pt", wgt);  
    histo.addVal(fabs(first.Eta() - second.Eta()), group,max, "DeltaEta", wgt); 
    histo.addVal(first.DeltaR(second),group,max, "DeltaR", wgt);  

    double dphiDijets = absnormPhi(second.Phi() - second.Phi());
    double dphi1 = normPhi(first.Phi() - theMETVector.Phi());
    double dphi2 = normPhi(second.Phi() - theMETVector.Phi());
    double alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histo.addVal(dphiDijets,group,max, "LeadSublDijetDphi", wgt); 
    //////histo.addVal(theMETVector.Pt(),dphiDijets, group,max, "");   //MetVsDiJetDeltaPhiLeadSubl
    /////histo.addVal(fabs(first.Eta() - second.Eta()), dphiDijets, group,max, "");   //DeltaEtaVsDeltaPhiLeadSubl

    histo.addVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), group,max, "R1", wgt);
    histo.addVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), group,max, "R2", wgt);
    histo.addVal(normPhi(first.Phi() - phiForMht), group,max, "Dphi1MHT", wgt); 
    histo.addVal(normPhi(second.Phi() - phiForMht), group,max, "Dphi2MHT", wgt);
    histo.addVal(dphi1,group,max, "Dphi1", wgt);
    histo.addVal(dphi2,group,max, "Dphi2", wgt);
    ///histo.addVal(dphi1,dphi2,group,max, "");   //Dphi1VsDphi2
    histo.addVal(alpha,group,max, "Alpha", wgt);


    //dijet info
  } else if(group == "FillDiJet") {
    double leaddijetmass = 0;
    double leaddijetpt = 0;
    double leaddijetdeltaR = 0;
    double leaddijetdeltaEta = 0;
    double etaproduct = -100;
    for(vec_iter it=goodParts[ival(CUTS::eDiJet)].begin(); it!=goodParts[ival(CUTS::eDiJet)].end(); it++) {
      int p1 = (*it) / _Jet->smearP.size();
      int p2 = (*it) % _Jet->smearP.size();
      TLorentzVector jet1 = _Jet->smearP.at(p1);
      TLorentzVector jet2 = _Jet->smearP.at(p2);
      TLorentzVector DiJet = jet1 + jet2;
	  
      if(DiJet.M() > leaddijetmass) {
	leaddijetmass = DiJet.M();
	etaproduct = (jet1.Eta() * jet2.Eta() > 0) ? 1 : -1;
      }
      leaddijetpt = (DiJet.Pt() > leaddijetpt) ? DiJet.Pt() : leaddijetpt;
      leaddijetdeltaEta = (fabs(jet1.Eta() - jet2.Eta()) > leaddijetdeltaEta) ? fabs(jet1.Eta() - jet2.Eta()) : leaddijetdeltaEta;
      leaddijetdeltaR = (jet1.DeltaR(jet2) > leaddijetdeltaR) ? jet1.DeltaR(jet2) : leaddijetdeltaR;

      histo.addVal(DiJet.M(), group,max, "Mass", wgt);
      histo.addVal(DiJet.Pt(), group,max, "Pt", wgt);
      histo.addVal(fabs(jet1.Eta() - jet2.Eta()), group,max, "DetlaEta", wgt);
      histo.addVal(absnormPhi(jet1.Phi() - jet2.Phi()), group,max, "DeltaPhi", wgt);
      histo.addVal(jet1.DeltaR(jet2), group,max, "DeltaR", wgt);
    }
  
    histo.addVal(leaddijetmass, group,max, "LeadMass", wgt);
    histo.addVal(leaddijetpt, group,max, "LeadPt", wgt);  
    histo.addVal(leaddijetdeltaEta, group,max, "LeadDeltaEta", wgt);
    histo.addVal(leaddijetdeltaR, group,max, "LeadDeltaR", wgt);
    histo.addVal(etaproduct, group,max, "LeadEtaProduct", wgt);


    ////diparticle stuff
  } else if(group == "FillDiMuon" || group == "FillDiTau" || group == "Muon1Tau1" || group == "Muon1Tau2" || group == "Muon2Tau1" || group == "Muon2Tau2") {  ///mumu/mutau/tautau
    Lepton* lep1 = NULL;
    Lepton* lep2 = NULL;
    CUTS ePos = fill_num[group];
    if(ePos == CUTS::eMuon1Tau1 || ePos == CUTS::eMuon1Tau2 || ePos == CUTS::eMuon2Tau1 || ePos == CUTS::eMuon2Tau2) {
      lep1 = _Muon; lep2 = _Tau;
    } else if(ePos == CUTS::eElec1Tau1 || ePos == CUTS::eElec1Tau2 || ePos == CUTS::eElec2Tau1 || ePos == CUTS::eElec2Tau2) {
      lep1 = _Electron; lep2 = _Tau;
    } else if(ePos == CUTS::eDiMuon) {
      lep1 = _Muon; lep2 = _Muon;
    } else if(ePos == CUTS::eDiTau) { lep1 = _Tau; lep2 = _Tau; 
    } else if (ePos == CUTS::eDiElec) { lep1 = _Electron; lep2 = _Electron; }

    
    for(vec_iter it=goodParts[ival(ePos)].begin(); it!=goodParts[ival(ePos)].end(); it++) {

      int p1= (*it) / _Gen->pt->size();
      int p2= (*it) % _Gen->pt->size();

      //    histo.addVal(Muon->smearP.at(mj).Pt(),Tau->smearP.at(tj).Pt());   //Muon1PtVsTau1Pt######################
      histo.addVal(lep1->smearP.at(p1).DeltaR(lep2->smearP.at(p2)), group,max, "DeltaR", wgt); 
      histo.addVal((lep1->smearP.at(p1).Pt() - lep1->smearP.at(p2).Pt()) / (lep1->smearP.at(p1).Pt() + lep2->smearP.at(p2).Pt()), group,max, "DeltaPtDivSumPt", wgt);   //Muon1Tau1DeltaPtDivSumPt
      histo.addVal(lep1->smearP.at(p1).Pt() - lep2->smearP.at(p2).Pt(), group,max, "DeltaPt", wgt);
      histo.addVal(cos(absnormPhi(lep2->smearP.at(p2).Phi() - lep1->smearP.at(p1).Phi())), group,max, "CosDphi", wgt);
      histo.addVal(absnormPhi(lep2->smearP.at(p2).Phi() - theMETVector.Phi()), group,max, "Part1MetDeltaPhi", wgt);
      ///      histo.addVal(absnormPhi(lep2->smearP.at(p2).Phi() - theMETVector.Phi()), cos(absnormPhi(lep2->smearP.at(p2).Phi() - lep1->smearP.at(p1).Phi())));   //Muon1MetDeltaPhiVsMuon1Tau1CosDphi
      histo.addVal(absnormPhi(lep1->smearP.at(p1).Phi() - theMETVector.Phi()), group,max, "Part2MetDeltaPhi", wgt);
      // if(CalculateTheDiTau4Momentum(lep1->smearP.at(p1),lep2->smearP.at(p2)).first) {
      // 	histo.addVal(DiParticleMass(lep1->smearP.at(p1),lep2->smearP.at(p2)), group,max, "ReconstructableMass");
      // } else {
      // 	histo.addVal(DiParticleMass(lep1->smearP.at(p1),lep2->smearP.at(p2)), group,max, "NotReconstructableMass");
      // }
      double PZeta = getPZeta(lep1->smearP.at(p1),lep2->smearP.at(p2));
      double PZetaVis = getPZetaVis(lep1->smearP.at(p1),lep2->smearP.at(p2));
      histo.addVal(calculateLeptonMetMt(lep2->smearP.at(p2)), group,max, "Part1MetMt", wgt);
      histo.addVal(calculateLeptonMetMt(lep1->smearP.at(p1)), group,max, "Part2MetMt", wgt); 
      histo.addVal(lep2->charge->at(p2) * lep1->charge->at(p1), group,max, "OSLS", wgt);  
      histo.addVal(PZeta, group,max, "PZeta", wgt); 
      histo.addVal(PZetaVis, group,max, "PZetaVis", wgt);  
      // histo.addVal(PZetaVis,PZeta, group,max, "Zeta2D");   //Muon1Tau1Zeta2D
      ///histo.addVal((_Muon1Tau1PZetaCutCoefficient * PZeta) + (_Muon1Tau1PZetaVisCutCoefficient * PZetaVis));   //Muon

      if ((goodParts[ival(CUTS::eR1stJet)].at(0) != -1) && (goodParts[ival(CUTS::eR2ndJet)].at(0) != -1)) {
	TLorentzVector TheLeadDiJetVect = _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].at(0)) + _Jet->smearP.at(goodParts[ival(CUTS::eR2ndJet)].size());
	histo.addVal(absnormPhi(lep1->smearP.at(p1).Phi() - TheLeadDiJetVect.Phi()), group,max, "Part1DiJetDeltaPhi", wgt);
	histo.addVal(absnormPhi(lep2->smearP.at(p2).Phi() - TheLeadDiJetVect.Phi()), group,max, "Part2DiJetDeltaPhi", wgt);
	////////#####histo.addVal(DiParticleMass(TheLeadDiJetVect, Muon->smearP.at(mj)+Tau->smearP.at(tj) ));   //Muon1Tau1DiJetReconstructableMass
      }
      // Fill(isZmm(Muon->smearP.at(mj)).first);   //Muon1Tau1_Muon1IsZmm
    }
  }
}


void Analyzer::initializePileupInfo(string MCHisto, string DataHisto) {
  // Filenames must be c_strings below. Here is the conversion from strings to c_strings
  // As you can see above cstr1 corresponds to MC and cstr2 corresponds to data.

  TFile *file1 = new TFile(MCHisto.c_str());
  TH1* histmc = static_cast<TH1*>(file1->Get("analyzeHiMassTau/NVertices_0"));
  if(!histmc) {throw std::runtime_error("failed to extract histogram");}
  for(int bin=0; bin<=(histmc->GetXaxis()->GetNbins() + 1); bin++) {
    hPUmc->SetBinContent(bin,histmc->GetBinContent(bin));
  }
  file1->Close();

  TFile* file2 = new TFile(DataHisto.c_str());
  TH1* histdata = static_cast<TH1*>(file2->Get("analyzeHiMassTau/NVertices_0"));
  if(!histdata) {throw std::runtime_error("failed to extract histogram");}
  for(int bin=0; bin<=(histdata->GetXaxis()->GetNbins() + 1); bin++) {
    hPUdata->SetBinContent(bin,histdata->GetBinContent(bin));
  }
  file2->Close();
}

double Analyzer::getPileupWeight(float ntruePUInt) {
  int bin;
  double MCintegral;
  double MCvalue;
  double Dataintegral;
  double Datavalue;

  // The probability that data (or MC) has N pileup interactions is value / integral
  // The ratio of the data and MC probability density functions gives us our pileup weights

  //std::cout << "Grabbing pileup info. " << std::endl;
  bin = hPUmc->GetBin(ntruePUInt+1);
  MCvalue = hPUmc->GetBinContent(bin);
  MCintegral = hPUmc->Integral();
  Datavalue = hPUdata->GetBinContent(bin);
  Dataintegral = hPUdata->Integral();

  // printouts for debugging
  //std::cout << "Number of true pileup interactions = " << ntruePUInt << std::endl;
  //std::cout << "Histogram bin, given the number of true pileup interactions = " << bin << std::endl;
  //std::cout << "MC PU probability density, given the number of true pileup interactions = " << MCvalue << std::endl;
  //std::cout << "Data PU probability density, given the number of true pileup interactions = " << Datavalue << std::endl;

  //std::cout << "Grabbing pileup weight. " << std::endl;
  //Ratio of normalized histograms in given bin

  return ((MCvalue * Dataintegral) != 0) ? (Datavalue * MCintegral) / (MCvalue * Dataintegral) : 1.0;
}





