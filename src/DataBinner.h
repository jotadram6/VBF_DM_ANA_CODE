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

class DataPiece {
 public:

  DataPiece() {};
  virtual ~DataPiece() {};
  virtual void write_histogram(vector<string>&, TFile*) {};
  virtual void bin(int, double, double) {};
  virtual void bin(int, double, double, double) {};
};


class Piece1D : public DataPiece {
 public:
  vector<float> data;
  double begin, end, width;
  int bins, Nfold;
  string name;
  TH1F histogram;

  Piece1D(string, int, double, double, int);
  ~Piece1D();
  void write_histogram(vector<string>&, TFile*);
  int get_bin(double);
  void bin(int, double, double);
};

class Piece2D : public DataPiece {
 public:
  vector<float> data;
  double beginx, endx, beginy, endy, widthx, widthy;
  int binx, biny, Nfold;
  string name;
  TH2F histogram;

  Piece2D(string, int, double, double, int, double, double, int);
  ~Piece2D();
    void write_histogram(vector<string>&, TFile*);
  int get_bin(double, bool);
  void bin(int, double, double, double);
};



class DataBinner {
public:
  DataBinner();
  ~DataBinner();
  void AddPoint(string,int, double, double);
  void AddPoint(string,int, double, double, double);
  void Add_Hist(string, string, int, double, double, int);
  void Add_Hist(string, string, int, double, double, int, double, double, int);
  void write_histogram(TFile*, vector<string>&);

private:
  unordered_map<string, DataPiece*> datamap;
  vector<string> order;
};

#endif
