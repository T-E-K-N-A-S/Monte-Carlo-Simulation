#include "file_stream.h"
#include "Headers.h"
#include "Parameters.h"
using namespace std;


void MatReader::read_phi(vector<vector < float> > &phi_mat,int mat_size,string file_name){
        vector< vector <float> > phi;
        vector <float> row;
        ifstream file;
        string line;
        string word;
        char delim = ',';
        file.open(file_name);
        if (file.is_open()){
            int i = 0;
            while(getline(file,line)){
                stringstream ss(line);                
                while(getline(ss,word,delim)){
                    row.push_back(stof(word));
                    i++;
                    if(i == mat_size){
                        phi.push_back(row);
                        row.clear();
                        i = 0;
                    }
                }
                
            }
        }
        else {
            cout << "unable to open file " + file_name << endl;
            exit(-1);
            
        }
        file.close();
        for(int i = 1 ; i < N-1; i++){
            for(int j = 1 ; j < N-1; j++){
                phi_mat[i][j] = phi[i-1][j-1];
            }
        }
}

void MatReader::read_spin(vector < vector <int> > &spin,int mat_size,string file_name){
        vector< vector <int> > occ;
        vector <int> row;
        
        ifstream file;
        string line;
        string word;
        char delim = ',';
        file.open(file_name);
        if (file.is_open()){
            int i = 0;
            while(getline(file,line)){
                stringstream ss(line);
                
                while(getline(ss,word,delim)){
                    if(word == "0"){
                        word = "-1";
                    }
                    row.push_back(stoi(word));
                    i++;
                    if(i == mat_size){
                        occ.push_back(row);
                        
                        row.clear();
                        i = 0;
                    }
                }
            }
        }
        else {
            cout << "unable to open file " + file_name << endl;
            exit(-1);
            
        }
        file.close();
        for(int i = 1 ; i < N-1; i++){
            for(int j = 1 ; j < N-1; j++){
                spin[i][j] = occ[i-1][j-1];
            }
        }
}


void MatWriter::write_spin(const vector <vector <int> > &spin, string file_name){
    ofstream f;
    f.open(file_name,ios::out);
    if (f.is_open()){
        for(auto i = 1 ; i != N-1 ; i++){
            for(auto j = 1 ; j != N-1 ; j++){
                f << spin[i][j] << ",";
            }
            f << std::endl;
        }

    }
    else {
        cout << "unable to open file : " +  file_name << endl;
        exit(-1);
    }
    f.close();

}


void MatWriter::write_phi(const vector < vector <float> > &phi,string file_name){
    ofstream f;
    f.open(file_name,ios::out);
    if (f.is_open()){
        for(auto i = 1 ; i != N-1 ; i++){
            for(auto j = 1 ; j != N-1 ; j++){
                f << phi[i][j] << ",";
            }
            f << std::endl;
        }        
    }
    else {
        cout << "unable to open file : " +  file_name << endl;
        exit(-1);
    }
    f.close();

}

