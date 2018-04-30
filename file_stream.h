#pragma once

#include "Headers.h"
using namespace std;
class MatReader{
public:
    static void read_phi(vector < vector <float> > & phi,int mat_size,string file_name);
        
    static void read_spin(vector < vector <int > > & spin,int mat_size,string file_name);
       

private:
    MatReader(){};
};


class MatWriter{
public:
    static void write_phi(const vector <vector < float> > &phi ,  string file_name);

    static void write_spin(const vector <vector < int > > &spin , string file_name); 

private:
    MatWriter(){};
};

