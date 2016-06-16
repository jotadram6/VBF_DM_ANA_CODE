#ifndef DATA_BINNER_H_
#define DATA_BINNER_H_

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

using namespace std;

struct DataPiece {
  vector<int> data;
  double begin;
  double width;
  int Nfold;
  
  DataPiece(double _begin, double _width, int _Nfold, int bins);
  int get_bin(double y);
  void bin(int folder, double y);
};



class DataBinner {
public:
  DataBinner();
  ~DataBinner();
  void AddPoint(string,int, double);
  void Add_Hist(string, int, double, double, int);

private:
  unordered_map<string, DataPiece*> datamap;
  vector<string> order;
};

#endif
