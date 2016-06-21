#ifndef DATA_BINNER_H_
#define DATA_BINNER_H_

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

using namespace std;

struct DataPiece {
  vector<int> data;
  double begin, end, width;
  int bins, Nfold;
  string name;
  TH1F histogram;

  DataPiece(string, int, double, double, int);
  ~DataPiece();
  void write_histogram(vector<string>&, TFile*);
  int get_bin(double);
  void bin(int, double);
};



class DataBinner {
public:
  DataBinner();
  ~DataBinner();
  void AddPoint(string,int, double);
  void Add_Hist(string, string, int, double, double, int);
  void write_histogram(TFile*, vector<string>&);

private:
  unordered_map<string, DataPiece*> datamap;
  vector<string> order;
};

#endif
