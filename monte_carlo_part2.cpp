#include "Headers.h"
#include "Parameters.h"
//#include "Clustering.h"
//#include "LatticeStuff.h"
#include "MonteCarlo.h"

using namespace std;

int monte_carlo_simulation(void)
{
	cout << "***** SIMULATED ANNEALING *****" << endl;
	cout << "Lattice size: " << LAT << endl;

	/* variables */
	int iseed,
		x1,
		ix,
		i_c,
		iy,
		j_c,
		imeas,
		iskip,
		nmcs,
		it,
		itemp,
		kconfig,
		k;
		

	vector <vector <int> > spin(N, vector <int>(N, 0));
	vector <vector <float> > phi(N, vector <float>(N, 0.0));
	vector <vector <float> > site_energy(N, vector <float>(N, 0.0));
	vector <float> temprature(TEMP, 0.0);
	vector <float> magnat(TEMP, 0.0);
	pair <int, int> rand_lattice_pt;

	
	float  phisum, boltz_dist, rand_var;
	float tot_en, tot_en_st, beta, mag, temp, delta_site_energy;

	// data file handling stuff 
	ofstream PhiMat, SpinMat, LogFile, f200, f300, f400, f12,file;
	PhiMat.open("phi_mat.txt");
	PhiMat << "i" << "\t" << "j" << "\t" << "phi[i][j]" << endl;
	SpinMat.open("spin_mat.txt");
	SpinMat << "i" << "\t" << "j" << "\t" << "spin[i][j]" << endl;
	LogFile.open("mylog_4.txt");
	file.open("magnat.csv");
	int st_flips = 0,pr_flips=0,rejec=0;
	iseed = 5;
	int mag2 = 0;

	/*	generate temperature	*/
	for (itemp = 0; itemp < TEMP; itemp++)
	{
		temprature[itemp] = 5- itemp*0.5;
		
	}


	//	dfh
	for (kconfig = 0; kconfig < NCONFIG; kconfig++)
	{
		//cout << "config " << kconfig << endl<< endl;

		generate_bimodal_dist(phi);	/*generates bimodal distribution over phi's*/
		save_matrix(PhiMat, phi);	/*save phi matrix using phimat ptr*/
		phisum = sum_of_matrix(phi);

		//cout << "creating spin mat " << endl;
		
		create_spin_mat(spin); /*make spin matrix*/
		//print_mat(spin, N);

		//	dfh for spin mat
		save_matrix(SpinMat, spin);
		
		//site_energy = siteEnergy(N, spin, phi, site_energy);
		//print_mat(phi, N);

		//cout << "creating hamiltonian" << endl;
		//	call hamiltonian()
		//hamiltonian(N, spin, phi, tot_en, tot_en_st);
		tot_en = cal_total_energy(spin, phi);
		//cout << "total energy\t" << tot_en << endl;

		LogFile << "initialen-toten" << "\t" << tot_en << "\t" << "encrystal=" << "\t" << 2 * N*N << "\t" << phisum << "\t" << "endisordercrystal=" << "\t" << -2 * N*N - phisum << "\t" << -2 * N*N + phisum << endl;
		LogFile << endl;
		LogFile << "kconfig" << "\t" << "temp" << "\t" << "\t" << "mag / NS" << "\t" << "toten" << endl;

		st_flips = 0; pr_flips = 0, rejec = 0;
		//cout << endl;
		for (itemp = 0; itemp < TEMP; itemp++)
		{
			//cout << "temp\t" << temprature[itemp] << endl;
			//cout << "No of Simul\t" << nmeass[itemp] << endl << endl;
			
			temp = temprature[itemp];
			beta = 1 / temp;
			//tot_en= 0;
			st_flips = 0; pr_flips = 0, rejec = 0;// , tot_en = 0;
			/***************************        Monte Carlo Simulation      *****************************/
			for (imeas = 0; imeas < MEAS; imeas++)

			{ 
				for (iskip = 0; iskip < SKIP; iskip++)

				{ // no of skip, sub epocs
				  
					for (nmcs = 0; nmcs < NS; nmcs++) // for n square, max rand guesses
					{
						
						rand_lattice_pt = generate_random_lattice_point();
						i_c = rand_lattice_pt.first;
						j_c = rand_lattice_pt.second;
						//cout << "i-> " << i_c << " j-> " << j_c << "\n";
		
						delta_site_energy =  -2 * energy_after_flip(spin,phi,i_c,j_c);
						if (delta_site_energy <= 0) // straight acceptance
						{
							//cout << "strt acceptance " << delta_site_energy << endl;
							tot_en += delta_site_energy;
							flip(spin, i_c, j_c);
							st_flips++;
							//cout<< "\t tot_en: " << tot_en << "," << cal_total_energy(spin, phi) << endl;
						}
						else
						{
							tot_en += metropolis(spin,delta_site_energy,beta,i_c,j_c,pr_flips,rejec);
							/* run metropolis algo */		
							//cout << "prob acceptance " << delta_site_energy << endl;
							//cout<< "\t tot_en: " << tot_en << "," << cal_total_energy(spin, phi) << endl;
						}
						//cout << "total energy\t" << tot_en << endl;

					}//cout << "total energy\t" << tot_en << endl;
					//system("pause");
				}
				//print_mat(spin,N);
				mag = 0;
				for (int i = 0; i < N; i++)
				{
					for (int j = 0; j < N; j++)
					{
						mag += spin[i][j];
					}
				}

				//cout << "mag_t" << mag << "\t" << endl;
				//cout << "del(mt)\t" << mag2 << "\t" << endl;
				mag2 = mag;

				f400.open("log_file" + to_string(kconfig) + ".txt", ofstream::app);
				f400 << "temp" << "\t" << "imeas" << "\t" << "mag / NS" << "\t" << "toten" << endl;
				f400 << temp << "\t" << imeas << "\t" << mag / NS << "\t\t" << tot_en << endl;
				f400 << endl;
				f400.close();

				/*cout << "\ntemp: "<<temp<<endl;*/
				cout << "#";
				//print_sign_mat(spin, N);
			}
			//cout << "\ntemp: " << temp << endl;
			//print_sign_mat(spin, N);
			//cluster call
			//clustering(spin);
			mag = 0;
			for (int i = 1; i < N-1; i++) {
				for (int j = 1; j < N-1; j++) {
					mag += spin[i][j];
				}
			}
			////cout << mag << endl;
			//LogFile << "kconfig" << "\t" << "temp" << "\t" << "\t" << "mag / NS" << "\t" << "toten" << endl;
			LogFile << kconfig << "\t" << temp << "\t" << "\t" << mag / NS << "\t" << tot_en << endl;
			cout	<< "\n\ncongif: "<<kconfig
					<< "\t temp: " << temp
					<< "\t mag: " <<mag
					<< "\t mag/site: "<< mag / pow(LAT,2)
					<< "\t tot_en: " << tot_en//<<","<<cal_total_energy(spin,phi)
					<< "\t max_en: "<<max_energy()
					<< "\n straight flips: "<<st_flips
					<< "\t prob flips: "<<pr_flips
					<< "\t rejections : " << rejec//float(st_flips)/float(pr_flips)
					<< "\t tot flips: " << st_flips + pr_flips + rejec
					<< endl;
			magnat.push_back(mag/(NS-2));
			file << abs(mag / pow((N - NDIM), 2)) << "\n";
		}
		st_flips = 0; pr_flips = 0,rejec=0;
		//int clusters = hoshen_kopelman(matrix, N,N);

		f200.open("spin_mat_2" + to_string(kconfig) + ".txt", ofstream::out);
		f200 << "i\tj\tspin[i][j]\n";

		f300.open("spin_mat_3" + to_string(kconfig) + ".txt", ofstream::out);
		f300 << "i\tj\tspin[i][j]\n";

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (spin[i][j] < 0) {
					f200 << i << "\t" << j << "\t" << spin[i][j] << endl;
				}
				else {

					f300 << i << "\t" << j << "\t" << spin[i][j] << endl;
				}
			}
		}
		f200.close();
		f300.close();

		f12 << kconfig << mag / NS << tot_en << endl;

		
	}
	f12.close();
	PhiMat.close();
	SpinMat.close();
	//printmat(phi);
	//cout << "st_flips\t" << flip << "\t" << endl;
	//printmat(site_en);
	//printmat1(spin);
	return 0;

}

int main(void) {

	monte_carlo_simulation();
	system("pause");
	return 0;
}