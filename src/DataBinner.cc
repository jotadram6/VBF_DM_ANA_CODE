#include "DataBinner.h"

using namespace std;

DataPiece::DataPiece(string name, int _bins, double _begin, double _end, int _Nfold) :  
  data(_Nfold*_bins, 0), begin(_begin), end(_end), bins(_bins), Nfold(_Nfold),
  histogram(name.c_str(), name.c_str(), _bins, _begin, _end) {
  
  width = (end - begin)/bins;
  
}

DataPiece::~DataPiece() {}
  

int DataPiece::get_bin(double y) {
  return (int)((y-this->begin)/this->width);
}

void DataPiece::bin(int folder, double y) {
  data[Nfold*folder + get_bin(y)]++;
}

void DataPiece::write_histogram(vector<string>& folders, TFile* outfile) {
  int entries;
  for(int i = 0; i<Nfold; i++) {
    
    outfile->cd(folders[i].c_str());
    entries = 0;
    for(int j = 0; j < bins; j++) {
      histogram.SetBinContent(j+1, data[i*Nfold + j]);
      entries += data[i*Nfold + j];
    }
    histogram.SetEntries(entries);
    histogram.Write();
  }
}


			      
 
DataBinner::DataBinner(){}

DataBinner::~DataBinner() {
  for(unordered_map<string, DataPiece*>::iterator it = datamap.begin(); it!=datamap.end(); it++) {
    delete it->second;
    it->second = NULL;
  }
}

void DataBinner::Add_Hist(string shortname, string fullname, int bin, double left, double right, int Nfolder) {
  datamap[shortname] = new DataPiece(fullname, bin, left, right, Nfolder);
  order.push_back(shortname);
}

void DataBinner::AddPoint(string name, int maxfolder, double value) {
  if(datamap[name] == NULL) return;
  for(int i=0; i < maxfolder; i++) {
      datamap[name]->bin(i,value);
  }
}

void DataBinner::write_histogram(TFile* outfile, vector<string>& folders) {
  for(vector<string>::iterator it = order.begin(); it != order.end(); it++) {
    datamap[*it]->write_histogram(folders, outfile);
  }
}


















