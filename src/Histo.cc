#include "Histo.h"

Histogramer::Histogramer(int _Npdf, string histname, string cutname, string outfilename): Npdf(_Npdf) {
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  read_cuts(cutname);
  NFolders = folders.size();

  read_hist(histname);
}

Histogramer::~Histogramer() {
  fill_histogram();

  for(vector<string>::iterator it = data_order.begin(); it != data_order.end(); it++) {
    delete data[*it];
    data[*it] = NULL;
  }
}


void Histogramer::read_hist(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");  

  if(!info_file) {
    cout << "ERROR: Didn't Read Histo File!" << endl;
    cout << filename << endl;
    exit(1);
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
	data[group] = new DataBinner();
	data_order.push_back(group); 
      }
    } else if(!accept) continue;
    else if(stemp.size() == 4) {
      string name = extractHistname(group, stemp[0]);
      data[group]->Add_Hist(name, stemp[0], stod(stemp[1]), stod(stemp[2]), stod(stemp[3]), NFolders);
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
  return histo;
}



void Histogramer::read_cuts(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");  

  if(!info_file) {
    cout << "ERROR: Didn't Read Histo File!" << endl;
    cout << filename << endl;
    exit(1);
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
	outfile->mkdir( name.c_str() );	
	folders.push_back(name);
	folder_num.push_back(i);
      } else if(stemp[1] == "0" && stemp[2] == "-1") continue;   ////remove unnecessary cuts
      cuts[name] = std::make_pair(stoi(stemp[1]),stoi(stemp[2]));
      cut_order.push_back(name);
      i++;
    }
  }

  info_file.close(); 
 }

void Histogramer::fill_histogram() {
  for(vector<string>::iterator it = data_order.begin(); it != data_order.end(); ++it) {
    data[*it]->write_histogram(outfile, folders);
  }
  outfile->Close();
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
    if(maxcut > folder_num[i]) maxFolder = folder_num[i];
    else break;
  }
  data[group]->AddPoint(histn, maxFolder, value);
}

// int main() {
//   Histogramer histo(1, "PartDet/Hist_entries.in", "PartDet/Cuts.in", "blah.root");
// }
    
