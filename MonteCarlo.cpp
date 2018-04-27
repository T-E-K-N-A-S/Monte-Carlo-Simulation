#include "Headers.h"	// >> contains all the req headerfiles
#include "Parameters.h" // >> contains all the req parameters
#include "MonteCarlo.h"


void print_sign_mat(std::vector <std::vector <int> >  & M, int len) {
	int i, j;
	for (i = 0; i < len; i++) {
		for (j = 0; j < len; j++) {
			if (M[i][j] == 0)std::cout << "" << "#";
			if (M[i][j] == 1)std::cout << "" << "+";
			if (M[i][j] == -1)std::cout << "" <<"-";
		}

		std::cout << "\n";
	}
}

/* Produces random floating-point values i, uniformly distributed on the interval [a, b),
that is, distributed according to the probability density function:
P(i|a,b) = 1 / ( b − a )   */
float uniformdist(void)
{
	std::random_device rd;
	std::mt19937 genarate(rd());
	std::uniform_real_distribution<> distribution(0, 1);

	return distribution(genarate);
}



void generate_bimodal_dist(std::vector <std::vector <float> > & phi)
{
	int i, j, n2, k;
	float var;
	for (i = 1; i < N - 1; i++)
	{
		for (j = 1; j < N - 1; j++)
		{
			var = uniformdist();
			////cout << aa << endl;
			if (var <= 0.5)
			{
				phi[i][j] = -delta;
			}
			else
			{
				phi[i][j] = delta;
			}
		}
	}
	//printmat(phi);
}

void create_spin_mat( std::vector<std::vector<int> >& spin)
{
	float a;


	for (int i = 1; i < N - 1; i++)
	{
		for (int j = 1; j < N - 1; j++)
		{
			a = 2 * uniformdist() - 1;
			if (a)
			{
				spin[i][j] = -1;
			}
			else
			{
				spin[i][j] = +1;
			}
		}
	}
	/*/// agumented matrix /// PBC
	for (int j = 1; j < N - 1; j++)
	{
	spin[N - 1][j] = spin[1][j];
	spin[0][j] = spin[N - 2][j];
	}
	for (int i = 1; i < N - 1; i++)
	{
	spin[i][N - 1] = spin[i][1];
	spin[i][0] = spin[i][N - 2];
	}*/

	/// agumented matrix /// NPBC
	for (int j = 1; j < N - 1; j++)
	{
		spin[N - 1][j] = 0;
		spin[0][j] = 0;
	}
	for (int i = 1; i < N - 1; i++)
	{
		spin[i][N - 1] = 0;
		spin[i][0] = 0;
	}

}

/* ------ H = sum_on(J.e_i.e_j.sig_i.sig_j) - B.sum_on(e_i.sig_i) --- eqn 6.4 pg 93---- */
/* ------- phi mat is B.e_i -------*/
void hamiltonian(std::vector< std::vector <int> > & spin, std::vector< std::vector <float> >& phi, float &tot_en, float &tot_en_st)
{
	float tot_en_2;
	tot_en = 0; 	tot_en_2 = 0;

	for (int i = 1; i < N - 1; i++)
	{
		for (int j = 1; j < N - 1; j++)
		{
			tot_en = tot_en - phi[i][j] * spin[i][j];

			/*			O		*/
			/*			|		*/
			/*	   O -- O -- O	*/
			/*			|		*/
			/*			O		*/

			tot_en_2 = tot_en_2 - spin[i][j] * (spin[i + 1][j] + spin[i - 1][j] + spin[i][j + 1] + spin[i][j - 1]);

		}
	}

	tot_en = tot_en + tot_en_2 / 2;
	tot_en_st = tot_en / (LAT*LAT);  /// error -solved ; old: tot_en/n*n
								 ////cout << "toten" << tot_en<< endl;


}

/* calculates total energy of lattice using hamilton's equation */
float cal_total_energy(std::vector<std::vector<int>>& spin, std::vector<std::vector<float>>& phi)
{
	float tot_en=0, partA=0,partB=0;
	

	for (int i = 1; i < N - 1; i++)
	{
		for (int j = 1; j < N - 1; j++)
		{
			/* partA accounts for disorder */
			partA = partA - phi[i][j] * spin[i][j]; 

			/*			O		*/
			/*			|		*/
			/*	   O -- O -- O   + disorder */
			/*			|		*/
			/*			O		*/

			/* partB accounts for nearest neibours */
			partB =  partB - spin[i][j] * (spin[i + 1][j] + spin[i - 1][j] + spin[i][j + 1] + spin[i][j - 1]);

		}
	}

	tot_en = partA + partB / 2;
	//tot_en_st = tot_en / (LAT*LAT);  /// error -solved ; old: tot_en/n*n
									 ////cout << "toten" << tot_en<< endl;

	return tot_en;
}

float energy_after_flip(std::vector<std::vector<int>>& spin, std::vector<std::vector<float>>& phi,int i, int j)
{
	float flip_en = 0;
	
	flip_en = -phi[i][j] * spin[i][j] /*disorder*/
			  -spin[i][j] * (spin[i + 1][j] + spin[i - 1][j] + spin[i][j + 1] + spin[i][j - 1]); /*neibours*/

	return flip_en;
}

void flip(std::vector<std::vector<int>>& spin, int i, int j)
{
	spin[i][j] *= -1;
}

float metropolis(std::vector<std::vector<int>>& spin,float del_en, float beta,int i,int j,int& pr_flips,int &rej)
{
	float rand_var = uniformdist();			// uniform
	float boltz = exp(-del_en*beta);		// boltzman
	float en = 0;
	if (rand_var <= boltz) // probabilistic acceptance
	{
		////cout << "probabilistic acceptance\n";
		//std::cout << "prob acceptance " << del_en << std::endl;
		en = del_en;
		flip(spin, i, j);
		//site_energy = update_HE(phi, i_c, j_c, site_energy, spin);	// update site energy
		pr_flips++;
		return en;
	}
	else
	{
		//std::cout << "rejection "  << std::endl;
		rej++;
		return 0;
	}
	
}

int max_energy()
{
	int max = LAT *LAT * 2;
	return max;
}
///  Stores site energy which corresponds to 4 neighbours in    site_en[][]
std::vector <std::vector <float> > siteEnergy(int n, std::vector <std::vector <int> >& spin, std::vector <std::vector <float> >& phi, std::vector <std::vector <float> >& site_en)
{
	for (int i = 1; i < n - 1; i++)
	{
		for (int j = 1; j < n - 1; j++)
		{	//new
			site_en[i][j] = -phi[i][j] * spin[i][j] - spin[i][j] * (spin[i + 1][j] + spin[i - 1][j] + spin[i][j - 1] + spin[i][j + 1]);
			//old
			//site_en[i][j] = - phi[i][j] - (spin[i+1][j] + spin[i-1][j] + spin[i][j-1] + spin[i][j+1]);

		}
	}

	return site_en;
}

// updates heartree energy; i dont know what the fuck heartree means, i am not sure its even spelled correctly.
std::vector <std::vector <float> > update_HE(std::vector <std::vector <float> >& phi, int i_c, int j_c, std::vector <std::vector <float> >& site_energy, std::vector <std::vector <int> >& spin)
{
	//new
	site_energy[i_c][j_c] = -phi[i_c][j_c] * spin[i_c][j_c] - spin[i_c][j_c] * (spin[i_c + 1][j_c] + spin[i_c - 1][j_c] + spin[i_c][j_c - 1] + spin[i_c][j_c + 1]);
	//	left
	if (i_c - 2 >= 0)
	{
		site_energy[i_c - 1][j_c] = -phi[i_c - 1][j_c] * spin[i_c - 1][j_c] - spin[i_c - 1][j_c] * (spin[i_c][j_c] + spin[i_c - 2][j_c] + spin[i_c - 1][j_c - 1] + spin[i_c - 1][j_c + 1]);
	}
	else
	{
		site_energy[i_c - 1][j_c] = -phi[i_c - 1][j_c] * spin[i_c - 1][j_c] - spin[i_c - 1][j_c] * (spin[i_c][j_c] + spin[i_c - 1][j_c - 1] + spin[i_c - 1][j_c + 1]);
	}
	//	right
	if (i_c + 2 < N)
	{
		site_energy[i_c + 1][j_c] = -phi[i_c + 1][j_c] * spin[i_c + 1][j_c] - spin[i_c + 1][j_c] * (spin[i_c + 2][j_c] + spin[i_c][j_c] + spin[i_c + 1][j_c - 1] + spin[i_c + 1][j_c + 1]);
	}
	else
	{
		site_energy[i_c + 1][j_c] = -phi[i_c + 1][j_c] * spin[i_c + 1][j_c] - spin[i_c + 1][j_c] * (spin[i_c][j_c] + spin[i_c + 1][j_c - 1] + spin[i_c + 1][j_c + 1]);
	}
	//	down
	if (j_c + 2 < N)
	{
		site_energy[i_c][j_c + 1] = -phi[i_c][j_c + 1] * spin[i_c][j_c + 1] - spin[i_c][j_c + 1] * (spin[i_c + 1][j_c + 1] + spin[i_c - 1][j_c + 1] + spin[i_c][j_c] + spin[i_c][j_c + 2]);
	}
	else
	{
		site_energy[i_c][j_c + 1] = -phi[i_c][j_c + 1] * spin[i_c][j_c + 1] - spin[i_c][j_c + 1] * (spin[i_c + 1][j_c + 1] + spin[i_c - 1][j_c + 1] + spin[i_c][j_c]);

	}
	//	up
	if (j_c - 2 >= 0)
	{
		site_energy[i_c][j_c - 1] = -phi[i_c][j_c - 1] * spin[i_c][j_c - 1] - spin[i_c][j_c - 1] * (spin[i_c + 1][j_c - 1] + spin[i_c - 1][j_c - 1] + spin[i_c][j_c - 2] + spin[i_c][j_c]);
	}
	else
	{
		site_energy[i_c][j_c - 1] = -phi[i_c][j_c - 1] * spin[i_c][j_c - 1] - spin[i_c][j_c - 1] * (spin[i_c + 1][j_c - 1] + spin[i_c - 1][j_c - 1] + spin[i_c][j_c]);
	}
	return site_energy;


}

std::pair <int, int> generate_random_lattice_point()
{
	int x = 0, y = 0;
	std::pair <int, int> latt_pt;
	x = int(uniformdist()*(LAT-1)+1);
	y = int(uniformdist()*LAT+1);
	//std::cout << x <<  "\t" <<y << std::endl;
	assert((x != LAT ) || (y != LAT+1));
	if ((x == LAT+1) || (y == LAT + 1)) latt_pt = generate_random_lattice_point();
	latt_pt.first = x;
	latt_pt.second = y;

	////cout << "rand_latt_pt: " << latt_pt.first << "," << latt_pt.second << endl;
	return latt_pt;
}


float sum_of_matrix(std::vector<std::vector<float> >& M)
{
	float sum = 0;
	for (int i = 0; i < M.size(); i++)
	{
		for (int j = 0; j < M[i].size(); j++)
		{
			sum += M[i][j]; // error: may be spin mat
		}
	}
	return sum;
}

int sum_of_matrix(std::vector<std::vector<int> >& M)
{
	float sum = 0;
	for (int i = 0; i < M.size(); i++)
	{
		for (int j = 0; j < M[i].size(); j++)
		{
			sum += M[i][j]; // error: may be spin mat
		}
	}
	return sum;
}



void print_mat(std::vector <std::vector <int> > &M, int length)
{
	int i, j;
	for (i = 0; i < length; i++)
	{
		for (j = 0; j < length; j++)
		{
			std::cout << M[i][j] << "\t";
		}
		std::cout << std::endl;
	}

}

void print_mat(std::vector <std::vector <float> > &M, int length)
{
	int i, j;
	for (i = 0; i < length; i++)
	{
		for (j = 0; j < length; j++)
		{
			std::cout << M[i][j] << "\t";
		}
		std::cout << std::endl;
	}

}

void save_matrix(std::ofstream &f, std::vector <std::vector <float> > &  M)
{
	for (auto i = 0; i < M.size(); i++)
	{
		for (auto j = 0; j < M[i].size(); j++)
		{
			f  << M[i][j] << ",";
		}
		f << std::endl;
	}
}

void save_matrix(std::ofstream &f, std::vector <std::vector <int> > &  M)
{
	for (auto i = 0; i < M.size(); i++)
	{
		for (auto j = 0; j < M[i].size(); j++)
		{
			f << M[i][j] << ",";
		}
		f << std::endl;
	}
}
