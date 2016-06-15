#include <vector.h>
#include 'Histo.h'

struct DataPiece {
  int* data;
  double begin;
  double width;
  int Nfold;

  int get_bin(double y) {
    return (int)((y-this.begin)/this.width);
  }

  void bin(int folder, double y) {
    data[Nfold*folder + get_bin(y)]++;
  }
  ~DataPiece() {
    delete[] data;
  }
};

DataBinner {
public:
  DataBinner(Histogramer);
  ~DataBinner();
  void AddPoint(int,int, double);

private:
  vector<DataPiece> dataArr;
  int arrsize;
}
