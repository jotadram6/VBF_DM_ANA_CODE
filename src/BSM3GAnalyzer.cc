void BSM3GAnalyzer::getInputs() {

  TEnv *params = new TEnv ("config_file");
  params->ReadFile ("BSM3GAnalyzer_CutParameters.in", kEnvChange);

  _TreatMuonsAsNeutrinos             = params->GetValue ("TreatMuonsAsNeutrinos", "string");


  _JetPtForMhtAndHt                  = params->GetValue ("JetPtForMhtAndHt", 1.0);
  _JetEtaForMhtAndHt                 = params->GetValue ("JetEtaForMhtAndHt", 1.0);
  _ApplyJetLooseIDforMhtAndHt        = params->GetValue ("ApplyJetLooseIDforMhtAndHt", "string");


  _Trigger1FirstRequirement             = params->GetValue ("Trigger1FirstRequirement", "string");
  _Trigger1SecondRequirement             = params->GetValue ("Trigger1SecondRequirement", "string");
  _Trigger2FirstRequirement             = params->GetValue ("Trigger2FirstRequirement", "string");
  _Trigger2SecondRequirement             = params->GetValue ("Trigger2SecondRequirement", "string");
  _CalculatePUSystematics            = params->GetValue ("CalculatePUSystematics", "string");


  _SmearTheJet                       = params->GetValue ("SmearTheJet", "string");
  _JetEnergyScaleOffset              = params->GetValue ("JetEnergyScaleOffset", 1.0);
  _MatchBToGen                       = params->GetValue ("MatchBToGen", "string");
  _DataHistos                        = params->GetValue ("DataHistos", "string");
  isData                             = params->GetValue ("isData", "string");
  _MCHistos                          = params->GetValue ("MCHistos", "string"); 


  string inputString;
  string inputType;
  


void BSM3GAnalyzer::analyze(TFile *theFile) {
  

  _totalEvents++;

  pdfWeightVector.clear();
  pdfWeightVector.push_back(1); // set the weight to 1 for the moment ...

  if((_CalculatePUSystematics == "1") && (isData == "0") ) { pu_weight = getPileupWeight(nTruePUInteractions); }
  else { pu_weight = 1.0; }

  deltaForMEx = 0;
  deltaForMEy = 0;


  if(_totalEvents == 1) {
    for(unsigned int i = 0 ; i < _TopologicalSelectionSequence.size(); i++) {
      string theDirectory = _TopologicalSelectionSequence[i];
      bookHistograms(theFile, theDirectory.c_str(), i);
    }
  }
  isrgluon_weight = pu_weight;
  
  getEventFlags();

  for(unsigned int i = 0 ; i < _EventSelectionSequence.size(); i++) {
    string theCut = _EventSelectionSequence[i];
    bool passedSelection = pass_i_EventSelectionSequence(i);
    for(unsigned int j = 0; j < _TopologicalSelectionSequence.size(); j++) {
      if(_TopologicalSelectionSequence[j] == theCut.c_str()) {
        for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++) {
          _hEvents[j][NpdfID]->Fill(0.0,isrgluon_weight * pdfWeightVector.at(NpdfID));
        }
        if(passedSelection == true) {
          //------Number of events passing cuts (numerator)
          for(unsigned int NpdfID = 0; NpdfID < pdfWeightVector.size();  NpdfID++) {
            _hEvents[j][NpdfID]->Fill(1.0,isrgluon_weight * pdfWeightVector.at(NpdfID));
          }
          fillHistograms(j);
        }
      }
    }
  }
  if(passEventSelectionSequence()) {
    _totalEventsPassingCuts++;
  }
}

void BSM3GAnalyzer::getEventFlags() {

  //---initiate event flags
  _EventFlag.clear();
  for (unsigned int i=0; i<_EventSelectionSequence.size(); i++) { _EventFlag.push_back(false); }

  //---Gen level tau requirements
  int test = 0;

  int nGenTaus = ExtractNumberOfGoodGen(15,2);
  if (nGenTaus>=_GenTauNmin) _EventFlag[_mapSelectionAlgoID["GenTauNmin"]] = true;
  if (nGenTaus<=_GenTauNmax) _EventFlag[_mapSelectionAlgoID["GenTauNmax"]] = true;

  //---Gen level top requirements
  int nGenTop = ExtractNumberOfGoodGen(6,2);
  if (nGenTop>=_GenTopNmin) _EventFlag[_mapSelectionAlgoID["GenTopNmin"]] = true;
  if (nGenTop<=_GenTopNmax) _EventFlag[_mapSelectionAlgoID["GenTopNmax"]] = true;

  //---Gen level electron requirements
  int nGenElectrons = ExtractNumberOfGoodGen(11,1);
  if (nGenElectrons>=_GenElectronNmin) _EventFlag[_mapSelectionAlgoID["GenElectronNmin"]] = true;
  if (nGenElectrons<=_GenElectronNmax) _EventFlag[_mapSelectionAlgoID["GenElectronNmax"]] = true;

  //---Gen level muon requirements
  int nGenMuons = ExtractNumberOfGoodGen(13,1);
  if (nGenMuons>=_GenMuonNmin) _EventFlag[_mapSelectionAlgoID["GenMuonNmin"]] = true;
  if (nGenMuons<=_GenMuonNmax) _EventFlag[_mapSelectionAlgoID["GenMuonNmax"]] = true;

  //---Gen level Z requirements
  int nGenZ = ExtractNumberOfGoodGen(23,2);
  if (nGenZ>=_GenZNmin) _EventFlag[_mapSelectionAlgoID["GenZNmin"]] = true;
  if (nGenZ<=_GenZNmax) _EventFlag[_mapSelectionAlgoID["GenZNmax"]] = true;

  //---Gen level W requirements
  int nGenW = ExtractNumberOfGoodGen(24,2);
  if (nGenW>=_GenWNmin) _EventFlag[_mapSelectionAlgoID["GenWNmin"]] = true;
  if (nGenW<=_GenWNmax) _EventFlag[_mapSelectionAlgoID["GenWNmax"]] = true;

  //---Gen level SM Higgs requirements
  int nGenSMHiggs = ExtractNumberOfGoodGen(25,2);
  if (nGenSMHiggs>=_GenSMHiggsNmin) _EventFlag[_mapSelectionAlgoID["GenSMHiggsNmin"]] = true;
  if (nGenSMHiggs<=_GenSMHiggsNmax) _EventFlag[_mapSelectionAlgoID["GenSMHiggsNmax"]] = true;

  //---Number of Good Vertices
  if (bestVertices>=_RecoVertexNmin) _EventFlag[_mapSelectionAlgoID["RecoVertexNmin"]] = true;
  if (bestVertices<=_RecoVertexNmax) _EventFlag[_mapSelectionAlgoID["RecoVertexNmax"]] = true;

  // ------Does the event pass the first trigger requirements?
  int nTriggersSatisfied = 0;
  if(passRecoTrigger1Cuts()) {nTriggersSatisfied++;}
  if (nTriggersSatisfied>=_RecoTriggers1Nmin) _EventFlag[_mapSelectionAlgoID["RecoTriggers1Nmin"]] = true;

  //---deltas for recalculation of MET (used when studying systematics or treating muons as neutrinos/taus)
  deltaForMEx = 0;
  deltaForMEy = 0;
  sumpxForMht = 0.0;
  sumpyForMht = 0.0;
  sumptForHt  = 0.0;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if( (smearedJetMomentumVector.at(j).Pt() > _JetPtForMhtAndHt) && (fabs(smearedJetMomentumVector.at(j).Eta()) < _JetEtaForMhtAndHt) ) {
      if(_ApplyJetLooseIDforMhtAndHt == "1") {
        if(passedLooseJetID(j)) {
          sumpxForMht = sumpxForMht - smearedJetMomentumVector.at(j).Px();
          sumpyForMht = sumpyForMht - smearedJetMomentumVector.at(j).Py();
          sumptForHt  = sumptForHt  + smearedJetMomentumVector.at(j).Pt();
        }
      } else {
        sumpxForMht = sumpxForMht - smearedJetMomentumVector.at(j).Px();
        sumpyForMht = sumpyForMht - smearedJetMomentumVector.at(j).Py();
        sumptForHt  = sumptForHt  + smearedJetMomentumVector.at(j).Pt();
      }
    }
  }
  if (sumpxForMht >= 0) {phiForMht = atan(sumpyForMht/sumpxForMht);}
  if (sumpxForMht < 0 && sumpyForMht >= 0) {phiForMht = atan(sumpyForMht/sumpxForMht) + TMath::Pi();}
  if (sumpxForMht < 0 && sumpyForMht < 0) {phiForMht = atan(sumpyForMht/sumpxForMht) - TMath::Pi();}

  //---Reco level muon1 requirements
  int nGoodCandidatesMuon1 = 0;
  for(int j = 0; j < Muon_pt->size(); j++) {
    if( ((passRecoMuon1Cuts(j)) || (passRecoMuon2Cuts(j))) && (_TreatMuonsAsNeutrinos == "1") ) {
      deltaForMEx += smearedMuonMomentumVector.at(j).Px();
      deltaForMEy += smearedMuonMomentumVector.at(j).Py();
    }
    if (!passRecoMuon1Cuts(j)) continue;
    nGoodCandidatesMuon1++;
  }
  if (nGoodCandidatesMuon1>=_RecoMuon1Nmin) _EventFlag[_mapSelectionAlgoID["RecoMuon1Nmin"]] = true;
  if (nGoodCandidatesMuon1<=_RecoMuon1Nmax) _EventFlag[_mapSelectionAlgoID["RecoMuon1Nmax"]] = true;

  //---Reco level muon2 requirements
  int nGoodCandidatesMuon2 = 0;
  for(int j = 0; j < Muon_pt->size(); j++) {
    if (!passRecoMuon2Cuts(j)) continue;
    nGoodCandidatesMuon2++;
  }
  if (nGoodCandidatesMuon2>=_RecoMuon2Nmin) _EventFlag[_mapSelectionAlgoID["RecoMuon2Nmin"]] = true;
  if (nGoodCandidatesMuon2<=_RecoMuon2Nmax) _EventFlag[_mapSelectionAlgoID["RecoMuon2Nmax"]] = true;

  //---Reco level electron1 requirements
  int nGoodCandidatesElectron1 = 0;
  for(int j = 0; j < patElectron_pt->size(); j++) {
    if (!passRecoElectron1Cuts(j)) continue;
    nGoodCandidatesElectron1++;
  }
  if (nGoodCandidatesElectron1>=_RecoElectron1Nmin) _EventFlag[_mapSelectionAlgoID["RecoElectron1Nmin"]] = true;
  if (nGoodCandidatesElectron1<=_RecoElectron1Nmax) _EventFlag[_mapSelectionAlgoID["RecoElectron1Nmax"]] = true;

  //---Reco level electron2 requirements
  int nGoodCandidatesElectron2 = 0;
  for(int j = 0; j < patElectron_pt->size(); j++) {
    if (!passRecoElectron2Cuts(j)) continue;
    nGoodCandidatesElectron2++;
  }
  if (nGoodCandidatesElectron2>=_RecoElectron2Nmin) _EventFlag[_mapSelectionAlgoID["RecoElectron2Nmin"]] = true;
  if (nGoodCandidatesElectron2<=_RecoElectron2Nmax) _EventFlag[_mapSelectionAlgoID["RecoElectron2Nmax"]] = true;

  //---Reco level tau1 requirements
  int nGoodCandidatesTau1 = 0;
  for(int j = 0; j < Tau_pt->size(); j++) {
    if (!passRecoTau1Cuts(j)) continue;
    nGoodCandidatesTau1++;
  }
  if (nGoodCandidatesTau1>=_RecoTau1Nmin) _EventFlag[_mapSelectionAlgoID["RecoTau1Nmin"]] = true;
  if (nGoodCandidatesTau1<=_RecoTau1Nmax) _EventFlag[_mapSelectionAlgoID["RecoTau1Nmax"]] = true;

  //---Reco level tau2 requirements
  int nGoodCandidatesTau2 = 0;
  for(int j = 0; j < Tau_pt->size(); j++) {
    if (!passRecoTau2Cuts(j)) continue;
    nGoodCandidatesTau2++;
  }
  if (nGoodCandidatesTau2>=_RecoTau2Nmin) _EventFlag[_mapSelectionAlgoID["RecoTau2Nmin"]] = true;
  if (nGoodCandidatesTau2<=_RecoTau2Nmax) _EventFlag[_mapSelectionAlgoID["RecoTau2Nmax"]] = true;

  //---recalculate MET
  double temppx = theMETVector.Px() + deltaForMEx;
  double temppy = theMETVector.Py() + deltaForMEy;
  double temppz = theMETVector.Pz();
  double temppt = TMath::Sqrt((temppx*temppx) + (temppy*temppy));
  TLorentzVector theTempMETVector;
  theTempMETVector.SetPxPyPzE(temppx,temppy,temppz,temppt);
  theMETVector = theTempMETVector;

  // ------Number of Good Jets
  int nGoodJets = 0;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if (!passRecoJet1Cuts(j)) continue;
    nGoodJets++;
  }
  if (nGoodJets>=_RecoJet1Nmin) _EventFlag[_mapSelectionAlgoID["RecoJet1Nmin"]] = true;
  if (nGoodJets<=_RecoJet1Nmax) _EventFlag[_mapSelectionAlgoID["RecoJet1Nmax"]] = true;

  nGoodJets = 0;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if (!passRecoJet2Cuts(j)) continue;
    nGoodJets++;
  }
  if (nGoodJets>=_RecoJet2Nmin) _EventFlag[_mapSelectionAlgoID["RecoJet2Nmin"]] = true;
  if (nGoodJets<=_RecoJet2Nmax) _EventFlag[_mapSelectionAlgoID["RecoJet2Nmax"]] = true;

  nGoodJets = 0;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if (!passRecoCentralJetCuts(j)) continue;
    nGoodJets++;
  }
  if (nGoodJets>=_RecoCentralJetNmin) _EventFlag[_mapSelectionAlgoID["RecoCentralJetNmin"]] = true;
  if (nGoodJets<=_RecoCentralJetNmax) _EventFlag[_mapSelectionAlgoID["RecoCentralJetNmax"]] = true;

  // ------Number of Good Leading Jets
  leadingjetpt = 0;
  leadingjeteta = -100;
  theLeadingJetIndex = -1;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if (!passRecoFirstLeadingJetCuts(j)) continue;
    if(smearedJetMomentumVector.at(j).Pt() > leadingjetpt) {
      leadingjetpt = smearedJetMomentumVector.at(j).Pt();
      leadingjeteta = smearedJetMomentumVector.at(j).Eta();
      theLeadingJetIndex = j;
    }
  }
  int nGoodFirstLeadingJets = 0;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if (j != theLeadingJetIndex) continue;
    nGoodFirstLeadingJets++;
  }
  if (nGoodFirstLeadingJets>=_RecoFirstLeadingJetNmin) _EventFlag[_mapSelectionAlgoID["RecoFirstLeadingJetNmin"]] = true;

  // ------Number of Good Second Leading Jets
  secondleadingjetpt = 0;
  secondleadingjeteta = -100;
  theSecondLeadingJetIndex = -1;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if (!passRecoSecondLeadingJetCuts(j)) continue;
    if((smearedJetMomentumVector.at(j).Pt() > secondleadingjetpt) && (j != theLeadingJetIndex)) {
      secondleadingjetpt = smearedJetMomentumVector.at(j).Pt();
      secondleadingjeteta = smearedJetMomentumVector.at(j).Eta();
      theSecondLeadingJetIndex = j;
    }
  }
  int nGoodSecondLeadingJets = 0;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if (j != theSecondLeadingJetIndex) continue;
    nGoodSecondLeadingJets++;
  }
  if (nGoodSecondLeadingJets>=_RecoSecondLeadingJetNmin) _EventFlag[_mapSelectionAlgoID["RecoSecondLeadingJetNmin"]] = true;

  // ------Number of Good b-Jets
  int nGoodBJets = 0;
  for(int j = 0; j < Jet_pt->size(); j++) {
    if (!passRecoBJetCuts(j)) continue;
    nGoodBJets++;
  }
  if (nGoodBJets>=_RecoBJetNmin) _EventFlag[_mapSelectionAlgoID["RecoBJetNmin"]] = true;
  if (nGoodBJets<=_RecoBJetNmax) _EventFlag[_mapSelectionAlgoID["RecoBJetNmax"]] = true;

  // ------Number of Good Susy Combinations (jet1+jet2+met combinations)
  int nGoodSusyCombinations = 0;
  if(passSusyTopologyCuts(theLeadingJetIndex,theSecondLeadingJetIndex)) {nGoodSusyCombinations++;}
  if (nGoodSusyCombinations>=_SusyCombinationsNmin) _EventFlag[_mapSelectionAlgoID["SusyCombinationsNmin"]] = true;

  // ------lepton+Met topology cuts
  nGoodCandidatesMuon1 = 0;
  nGoodCandidatesElectron1 = 0;
  nGoodCandidatesTau1 = 0;
  nGoodCandidatesMuon2 = 0;
  nGoodCandidatesElectron2 = 0;
  nGoodCandidatesTau2 = 0;
  for(int j = 0; j < Muon_pt->size(); j++) {
    if (!passRecoMuon1Cuts(j)) continue;
    if (!passRecoMuon1MetTopologyCuts(j)) continue;
    nGoodCandidatesMuon1++;
  }
  for(int j = 0; j < Muon_pt->size(); j++) {
    if (!passRecoMuon2Cuts(j)) continue;
    if (!passRecoMuon2MetTopologyCuts(j)) continue;
    nGoodCandidatesMuon2++;
  }
  for(int j = 0; j < patElectron_pt->size(); j++) {
    if (!passRecoElectron1Cuts(j)) continue;
    if (!passRecoElectron1MetTopologyCuts(j)) continue;
    nGoodCandidatesElectron1++;
  }
  for(int j = 0; j < patElectron_pt->size(); j++) {
    if (!passRecoElectron2Cuts(j)) continue;
    if (!passRecoElectron2MetTopologyCuts(j)) continue;
    nGoodCandidatesElectron2++;
  }
  for(int j = 0; j < Tau_pt->size(); j++) {
    if (!passRecoTau1Cuts(j)) continue;
    if (!passRecoTau1MetTopologyCuts(j)) continue;
    nGoodCandidatesTau1++;
  }
  for(int j = 0; j < Tau_pt->size(); j++) {
    if (!passRecoTau2Cuts(j)) continue;
    if (!passRecoTau2MetTopologyCuts(j)) continue;
    nGoodCandidatesTau2++;
  }
  if (nGoodCandidatesMuon1>=_RecoMuon1MetTopologyNmin) _EventFlag[_mapSelectionAlgoID["RecoMuon1MetTopologyNmin"]] = true;
  if (nGoodCandidatesMuon1<=_RecoMuon1MetTopologyNmax) _EventFlag[_mapSelectionAlgoID["RecoMuon1MetTopologyNmax"]] = true;
  if (nGoodCandidatesMuon2>=_RecoMuon2MetTopologyNmin) _EventFlag[_mapSelectionAlgoID["RecoMuon2MetTopologyNmin"]] = true;
  if (nGoodCandidatesMuon2<=_RecoMuon2MetTopologyNmax) _EventFlag[_mapSelectionAlgoID["RecoMuon2MetTopologyNmax"]] = true;
  if (nGoodCandidatesElectron1>=_RecoElectron1MetTopologyNmin) _EventFlag[_mapSelectionAlgoID["RecoElectron1MetTopologyNmin"]] = true;
  if (nGoodCandidatesElectron1<=_RecoElectron1MetTopologyNmax) _EventFlag[_mapSelectionAlgoID["RecoElectron1MetTopologyNmax"]] = true;
  if (nGoodCandidatesElectron2>=_RecoElectron2MetTopologyNmin) _EventFlag[_mapSelectionAlgoID["RecoElectron2MetTopologyNmin"]] = true;
  if (nGoodCandidatesElectron2<=_RecoElectron2MetTopologyNmax) _EventFlag[_mapSelectionAlgoID["RecoElectron2MetTopologyNmax"]] = true;
  if (nGoodCandidatesTau1>=_RecoTau1MetTopologyNmin) _EventFlag[_mapSelectionAlgoID["RecoTau1MetTopologyNmin"]] = true;
  if (nGoodCandidatesTau1<=_RecoTau1MetTopologyNmax) _EventFlag[_mapSelectionAlgoID["RecoTau1MetTopologyNmax"]] = true;
  if (nGoodCandidatesTau2>=_RecoTau2MetTopologyNmin) _EventFlag[_mapSelectionAlgoID["RecoTau2MetTopologyNmin"]] = true;
  if (nGoodCandidatesTau2<=_RecoTau2MetTopologyNmax) _EventFlag[_mapSelectionAlgoID["RecoTau2MetTopologyNmax"]] = true;

// ------------ method called once each job just after ending the event loop  ------------
void BSM3GAnalyzer::endJob(TFile *theOutFile) {

  printEfficiency();  

  theOutFile->cd();
  for(unsigned int i = 0 ; i < _TopologicalSelectionSequence.size(); i++) {
    string theDirectory = _TopologicalSelectionSequence[i];
    writeHistograms(theOutFile, theDirectory.c_str(), i);
  }
  theOutFile->Close();

}





bool BSM3GAnalyzer::passDiJetTopologyCuts(int nobj1, int nobj2) {
  // ----Separation cut between jets (remove overlaps)
  if (_DoDiJetDiscrByDeltaR == "1") {if(smearedJetMomentumVector.at(nobj1).DeltaR(smearedJetMomentumVector.at(nobj2)) < _DiJetDeltaRCut) {return false;}}
  if (_DoDiJetDiscrByDeltaEta == "1") {
    if(fabs(smearedJetMomentumVector.at(nobj1).Eta() - smearedJetMomentumVector.at(nobj2).Eta()) < _DiJetMinDeltaEtaCut) {return false;}
    if(fabs(smearedJetMomentumVector.at(nobj1).Eta() - smearedJetMomentumVector.at(nobj2).Eta()) > _DiJetMaxDeltaEtaCut) {return false;}
  }
  if (_DoDiJetDiscrByDeltaPhi == "1") {
    if(abs(normalizedPhi(smearedJetMomentumVector.at(nobj1).Phi() - smearedJetMomentumVector.at(nobj2).Phi())) < _DiJetMinDeltaPhiCut) {return false;}
    if(abs(normalizedPhi(smearedJetMomentumVector.at(nobj1).Phi() - smearedJetMomentumVector.at(nobj2).Phi())) > _DiJetMaxDeltaPhiCut) {return false;}
  } 
  if (_DoDiJetDiscrByOSEta == "1") {
    if((smearedJetMomentumVector.at(nobj1).Eta() * smearedJetMomentumVector.at(nobj2).Eta()) >= 0) {return false;}
  }
  // ----Require both legs to be almost back-to-back in phi
  if (_DoDiJetDiscrByCosDphi == "1") {
    if(cos(abs(normalizedPhi(smearedJetMomentumVector.at(nobj1).Phi() - smearedJetMomentumVector.at(nobj2).Phi()))) > _DiJetCosDphiMaxCut) {return false;}
    if(cos(abs(normalizedPhi(smearedJetMomentumVector.at(nobj1).Phi() - smearedJetMomentumVector.at(nobj2).Phi()))) < _DiJetCosDphiMinCut) {return false;}
  }
  // ----Mass window requirement
  if (_DoDiscrByDiJetMassReco == "1") {
    if( (CalculateTheDiJet4Momentum(smearedJetMomentumVector.at(nobj1),smearedJetMomentumVector.at(nobj2)).second.M() < _DiJetMassMinCut) || (CalculateTheDiJet4Momentum(smearedJetMomentumVector.at(nobj1),smearedJetMomentumVector.at(nobj2)).second.M() > _DiJetMassMaxCut) ) {return false;}
  }
  return true;
}

