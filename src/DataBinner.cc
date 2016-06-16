#include "DataBinner.h"

using namespace std;

DataPiece::DataPiece(double _begin, double _width, int _Nfold, int bins) 
  : begin(_begin), width(_width), Nfold(_Nfold) {

  data.resize(Nfold*bins);
}

int DataPiece::get_bin(double y) {
  return (int)((y-this->begin)/this->width);
}

void DataPiece::bin(int folder, double y) {
  data[Nfold*folder + get_bin(y)]++;
}


DataBinner::DataBinner(){}

DataBinner::~DataBinner() {
  for(unordered_map<string, DataPiece*>::iterator it = datamap.begin(); it!=datamap.end(); it++) {
    delete it->second;
    it->second = NULL;
  }
}

void DataBinner::Add_Hist(string name, int bin, double left, double right, int Nfolder) {
  double width = (right - left) / bin;
  DataPiece* temp = new DataPiece(left, width, Nfolder, bin);
  datamap[name] = temp;
  order.push_back(name);
}

void DataBinner::AddPoint(string name, int maxfolder, double value) {
  if(datamap[name] == NULL) return;
  for(int i=0; i < maxfolder; i++) {
      datamap[name]->bin(i,value);
  }
}

