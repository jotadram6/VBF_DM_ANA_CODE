#include "Histo.h"

typedef std::pair<string, std::array<double, 3> > info_el;
typedef std::vector<std::pair<string, std::array<double, 3> > > info_vec;
typedef info_vec::iterator info_it;


Histogramer::Histogramer(int _Npdf, string fname): Npdf(_Npdf) {
  outfile = new TFile(fname.c_str(), "RECREATE");

  read_hist();
  read_cuts();

  NFolders = folders.size();

  // for(vector<string>::iterator it=folders.begin(); it !=folders.end(); it++) {
  //   std::cout << *it << std::endl;
  // }
  // for(info_it iter=Generator_info.begin(); iter != Generator_info.end(); iter++) {
  //   std::cout << (*iter).first << " " << (*iter).second[0] << " " << (*iter).second[1] << " " << (*iter).second[2] << std::endl;
  // }
  
  for(int i = 0; i < NFolders; i++) {
    string directory = folders[i];
    //TDirectory *subDir = 
    outfile->mkdir( directory.c_str() );
  }
  for(info_it iter=Generator_info.begin(); iter != Generator_info.end(); iter++) {
    std::vector<TH1*> tmpVec;
    Generator_Histogram[(*iter).first] = tmpVec;

    for(int i = 0; i < NFolders; i++) {
      outfile->cd(folders[i].c_str() );

      for(int j = 0; j < Npdf; j++) {
	string name = (*iter).first + "_" + to_string(j);
	TH1F* hist = new TH1F(name.c_str(), name.c_str(), (*iter).second[0], (*iter).second[1], (*iter).second[2]);
	Generator_Histogram[(*iter).first].push_back(hist);
      }
    }
  }


}

Histogramer::~Histogramer() {
  write_histogram();
  outfile->Close();
}


void Histogramer::read_hist() {
  ifstream info_file("Hist_entries.in");
  
  if(!info_file) return;

  string line;
  string current_area = "";
  bool accept = false;
  while(getline(info_file, line)) {
    stringstream ss(line);
    double bin, left, right, extra;
    string name;
    if(!(ss >> name)) {
      continue;
    }
    if(!(ss >> bin >> left)) {
      current_area = name;
      accept = int(bin);
    } else if(!(ss >> right >> extra) && accept) {
      std::array<double, 3> tmpVals { bin, left, right };
      info_el tmpPair(name, tmpVals);
      Generator_info.push_back(tmpPair);
     }
  }
}

void Histogramer::read_cuts() {
  ifstream info_file("Cuts.in");
  
  if(!info_file) return;

  string line;
  while(getline(info_file, line)) {
    stringstream ss(line);
    double min, max;
    string name;
    if(!(ss >> name)) continue;
    if(name[0] == '/' && name[1] == '/') continue;

    if((ss >> min >> max)) {
      if(name[0]=='*' && name[1]=='*' && name[2]=='*') {
	name.erase(0,3);
	folders.push_back(name);
      }
      cuts.push_back(name);
      cut_range.push_back(std::make_pair(min,max));
    }
  }
}

int Histogramer::index(int x, int y) {
  return x*NFolders + y;
}

void Histogramer::write_histogram() {
  for(info_it iter=Generator_info.begin(); iter != Generator_info.end(); iter++) {
    for(int i = 0; i < NFolders; i++) {
      outfile->cd( folders[i].c_str() );
      
      for(int j = 0; j < Npdf; j++) {
	Generator_Histogram[(*iter).first][index(i,j)]->Write();
      }
    }
  }  
}

int Histogramer::get_Nhists(){
  return Generator_info.size();
}

double Histogramer::get_start(int hist) {
  return Generator_info[hist].second[1];
}

double Histogramer::get_width(int hist) {
  return (Generator_info[hist].second[2] - Generator_info[hist].second[1]) / Generator_info[hist].second[0];
}

int Histogramer::get_nbins(int hist) {
    return Generator_info[hist].second[0];
}


int main() {
  Histogramer test(4,"blah.root");
  std::cout << test.get_Nhists() << std::endl;
  int num = test.get_Nhists();
  for( int i=0; i<num; i++) {
    std::cout << "start " << test.get_start(i) << std::endl;
    std::cout << "width " << test.get_width(i) << std::endl;
    std::cout << "nbins " << test.get_nbins(i) << std::endl;
  }
  
}
