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
#include <string>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <regex>
#include "DataBinner.h"
#include "tokenizer.hpp"


using namespace std;

class Histogramer {

 public:
  Histogramer(int, string, string, string);
  ~Histogramer();
  
  vector<int> get_folders();
  
  unordered_map<string,pair<int,int>>* get_cuts();
  vector<string>* get_order();
  vector<string>* get_groups();
  void addVal(double, string, int, string, double);

 private:
  TFile * outfile;

  int NFolders;
  int Npdf;

  unordered_map<string, pair<int,int>> cuts;
  vector<string> cut_order;
  vector<string> folders;
  vector<int> folder_num;
  unordered_map<string, DataBinner*> data;
  vector<string> data_order;

  void read_hist(string);
  void read_cuts(string);
  void fill_histogram();
  
  string extractHistname(string, string);
};

#endif

