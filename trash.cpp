#include "file_stream.h"
//#include "Headers.h"
using namespace std;

int main(void){
    string csv_path = "csvs";
    MatReader::read_phi(128,csv_path + "Si-1-13.csv");
    return 0;
}