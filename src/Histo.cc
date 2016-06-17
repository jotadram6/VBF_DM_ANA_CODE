#include "Histo.h"

typedef std::pair<string, std::array<double, 3> > info_el;
typedef std::vector<std::pair<string, std::array<double, 3> > > info_vec;
typedef info_vec::iterator info_it;


Histogramer::Histogramer(int _Npdf, string histname, string cutname, string outfilename): Npdf(_Npdf) {
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  read_cuts(cutname);
  NFolders = folders.size();

  read_hist(histname);
  
  // for(vector<string>::iterator it=folders.begin(); it !=folders.end(); it++) {
  //   std::cout << *it << std::endl;
  // }
  // for(info_it iter=Generator_info.begin(); iter != Generator_info.end(); iter++) {
  //   std::cout << (*iter).first << " " << (*iter).second[0] << " " << (*iter).second[1] << " " << (*iter).second[2] << std::endl;
  // }

  for(int i = 0; i < NFolders; i++) {
    string directory = folders[i].first;
    outfile->mkdir( directory.c_str() );
  }
  for(info_it iter=Generator_info.begin(); iter != Generator_info.end(); iter++) {
    std::vector<TH1*> tmpVec;
    Generator_Histogram[(*iter).first] = tmpVec;
    for(int i = 0; i < NFolders; i++) {
      outfile->cd(folders[i].first.c_str() );

      for(int j = 0; j < Npdf; j++) {
  	string name = (*iter).first + "_" + to_string(j);
  	TH1F* hist = new TH1F(name.c_str(), name.c_str(), (*iter).second[0], (*iter).second[1], (*iter).second[2]);
  	Generator_Histogram[(*iter).first].push_back(hist);
      }
    }
  }


}

Histogramer::~Histogramer() {
  for(vector<string>::iterator it = data_order.begin(); it != data_order.end(); it++) {
    delete data[*it];
    data[*it] = NULL;
  }
  write_histogram();
  outfile->Close();
}


void Histogramer::read_hist(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");  

  if(!info_file) {
    cout << "no read histo!" << endl;
    return;
  }

  vector<string> stemp;
  string group,line;
  bool accept = false;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }

    if(stemp.size() == 0) continue;
    else if(stemp.size() == 2) {
      group = stemp[0];
      accept = stoi(stemp[1]);
      if(accept) {
	DataBinner* dbtemp = new DataBinner();
	data[group] = dbtemp;
	data_order.push_back(group); 
	cout << endl << group << endl;
      }
    } else if(!accept) continue;
    else if(stemp.size() == 4) {
      std::array<double, 3> tmpVals {stod(stemp[1]), stod(stemp[2]), stod(stemp[3])};
      info_el tmpPair(stemp[0], tmpVals);
      Generator_info.push_back(tmpPair);
      string name = extractHistname(group, stemp[0]);
      data[group]->Add_Hist(name,stod(stemp[1]), stod(stemp[2]), stod(stemp[3]), NFolders);
    }
  }

  info_file.close(); 
}

string Histogramer::extractHistname(string group, string histo) {
  regex reg ("((Tau|Muon|Electron)+(1|2)+)");
  smatch m;

  string stringkey = group.erase(0,4);
  regex first (stringkey+"(_)?");
  histo = regex_replace(histo,first, "");
  if(stringkey.find("Di") != string::npos) {
    stringkey=stringkey+"1"+stringkey+"2";
  }

  int i=1;
  while(regex_search(stringkey,m,reg)) {
    regex key(m[0].str());
    histo = regex_replace(histo,key,"Part"+to_string(i));
    stringkey = m.suffix().str();
    i++;
  }
  cout << histo << endl;
  return histo;
}
  



void Histogramer::read_cuts(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");  

  if(!info_file) {
    cout << "no read histo!" << endl;
    return;
  }

  vector<string> stemp;
  string name,line;
  int i = 0;

  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }

    if(stemp.size() == 3) {
      name = stemp[0];
      if(name[0]=='*' && name[1]=='*' && name[2]=='*') {
	name.erase(0,3);
	folders.push_back(make_pair(name,i));
      }
      cuts[name] = std::make_pair(stoi(stemp[1]),stoi(stemp[2]));
      cut_order.push_back(name);
      i++;
    }
  }

  info_file.close(); 
 }

int Histogramer::index(int x, int y) {
  return x*NFolders + y;
}

void Histogramer::write_histogram() {
  for(info_it iter=Generator_info.begin(); iter != Generator_info.end(); iter++) {
    for(int i = 0; i < NFolders; i++) {
      outfile->cd( folders[i].first.c_str() );
      
      for(int j = 0; j < Npdf; j++) {
	Generator_Histogram[(*iter).first][index(i,j)]->Write();
      }
    }
  }  
}





unordered_map<string,pair<int,int>>* Histogramer::get_cuts() {
  return &cuts;
}

vector<string>* Histogramer::get_order() {
  return &cut_order;
}

vector<string>* Histogramer::get_groups() {
  return &data_order;
}


void Histogramer::addVal(double value, string group, int maxcut, string histn) {
  int maxFolder=0;
  for(int i = 0; i < NFolders; i++) {
    if(maxFolder > folders[i].second) maxFolder = folders[i].second;
    else break;
  }
  data[group]->AddPoint(histn, maxFolder, value);
}

// int main() {
//   Histogramer* test = new Histogramer(4,"PartDet/Hist_entries.in", "PartDet/Cuts.in","blah.root");
//   delete test;
  
// }
