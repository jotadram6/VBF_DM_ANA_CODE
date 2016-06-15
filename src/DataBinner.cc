
DataBinner::DataBinner(Histogramer histo) {
  arrsize = histo.get_Nhists();
  
  for(int i=0; i < arrsize; i++) {
    DataPiece temp = {};
    temp.begin = histo.get_start(i);
    temp.width = histo.get_width(i);
    data = new int[histo.Nfolders * histo.get_nbins(int)];
    dataArr.push_back(temp);
  }
}

void AddPoint(int histogram, int folder, double value) {
  dataArr[histogram].bin(folder,value);
}


