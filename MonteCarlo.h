#pragma once

#include "Headers.h"
/*	probability */
float uniformdist(void);
void generate_bimodal_dist(std::vector <std::vector <float> >& phi);
std::pair <int, int> generate_random_lattice_point();

/*	core */
void create_spin_mat(std::vector<std::vector<int> >& spin);
void hamiltonian(std::vector< std::vector <int> >& spin, std::vector< std::vector <float> >& phi, float &tot_en, float &tot_en_st);
std::vector <std::vector <float> > siteEnergy(int n, std::vector <std::vector <int> >& spin, std::vector <std::vector <float> >& phi, std::vector <std::vector <float> >& site_en);
std::vector <std::vector <float> > update_HE(std::vector <std::vector <float> >& phi, int i_c, int j_c, std::vector <std::vector <float> >& site_energy, std::vector <std::vector <int> >& spin);
float sum_of_matrix(std::vector <std::vector <float> > &M);
int sum_of_matrix(std::vector<std::vector<int> >& M);
float cal_total_energy(std::vector< std::vector <int> > & spin, std::vector< std::vector <float> >& phi);
float energy_after_flip(std::vector< std::vector <int> > & spin, std::vector< std::vector <float> >& phi,int i,int j);
void flip(std::vector< std::vector <int> > & spin, int i, int j);
float metropolis(std::vector< std::vector <int> > & spin, float del_en, float beta, int i, int j,int& pr_flips,int &rej);
int max_energy();
/*	overheads	*/

void save_matrix(std::ofstream &f, std::vector <std::vector <int> > &  M);
void save_matrix(std::ofstream &f, std::vector <std::vector <float> > &  M);

void print_mat(std::vector <std::vector <int> > & M, int length);
void print_mat(std::vector <std::vector <float> > & M, int length);
void print_sign_mat(std::vector <std::vector <int> >  & M, int len);

