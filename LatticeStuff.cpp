#include "Headers.h"
#include "Parameters.h" // >> contains all the req parameters
 
void create_Bmat_bimodal(std::vector <float>& Bmat,float del)
{
	long i = 0, var = 0, count1 = 0, count2 = 0;
	std::uniform_int_distribution<int>  distr(0, 10000);
	std::random_device  rand_dev;
	std::mt19937        generator(rand_dev());
	std::cout << "B mtrix ->" << tab;
	for (i = 0; i < N; i++)
	{
		var = distr(generator);
		if (var < 5000)
		{
			Bmat[i] = del;
			count1++;
		}
		else
		{
			Bmat[i] = -del;
			count2++;
		}
		//std::cout << Bmat[i] << d;
	}	std::cout << std::endl << "counts : " << count1 << d << count2 << std::endl;


}

//
void update_Bmat_bimodal(std::vector <float>& Bmat, float del)
{
	for (int i = 0; i < N; i++)
	{
		
		if (Bmat[i] > 0)
		{
			Bmat[i] = del;
			
		}
		else
		{
			Bmat[i] = -del;
			
		}
		//std::cout << Bmat[i] << d;
	}
}
// http://www.cplusplus.com/reference/random/normal_distribution/
void create_Bmat_gaussian(std::vector <float>& Bmat,float del)
{
	double var = 0;
	long i = 0,  count1 = 0, count2 = 0;
	std::normal_distribution<float>  distr(0, float(del)/10);
	std::random_device  rand_dev;
	std::mt19937        generator(rand_dev());
	std::cout << "B mtrix ->" << tab;
	for (i = 0; i < N; i++)
	{
		//Bmat[i] = distr(generator);// floor(rand() / 1000);
		var = (distr(generator));
		Bmat[i] = var;
		if (var < 0)
		{
			//Bmat[i] = w;
			count1++;
		}
		else
		{
			//Bmat[i] = -w;
			count2++;
		}
		//std::cout << Bmat[i] << d;
	}	std::cout << std::endl << "counts : " << count1 << d << count2 << std::endl;


}

void create_Exmat(std::vector < std::vector <int> >& Exmat, std::vector <int> & latt)
{
	int i = 0;
	for (i = 0; i < N; i++)
	{
		if (((i + 1) % VER != 0))
		{
			Exmat[i][i + 1] = 1 * latt[i] * latt[i + 1];
			//cout << Exmat[i][i + 1];
		}
		if (i + VER < N)
		{
			Exmat[i][i + VER] = 1 * latt[i] * latt[i + VER];
		}
		//else 
	}
}


void create_CapacityMat(std::vector < std::vector <float> >& CapacityMat, std::vector < std::vector <int> >& Exmat)
{
	for (int i = 1; i <= N; i++)
	{
		//cout << "\n";
		for (int j = 1; j <= N; j++)
		{
			if (i < j)
			{
				CapacityMat[i][j] = 4 * Exmat[i - 1][j - 1] ;

			}
			else
				CapacityMat[i][j] = 0;
			//cout << CapacityMat[i][j] << " ";
		}
	}
}


void create_Wmat(std::vector < float >& Wmat,std::vector < std::vector <float> >& CapacityMat, std::vector < float >& Bmat, float del )
{
	int i = 0, j = 0, cap = 0;
	for (i = 1; i <= N; i++)
	{
		for (j = 1; j <= N; j++)
		{
			cap += CapacityMat[i][j] - CapacityMat[j][i];
			//cout << cap<<"\t"; 
		}
		Wmat[i - 1] = -2 * Bmat[i - 1] * del - cap / 2;
		cap = 0;
		//std::cout <<Wmat[i-1]<< std::endl;
	}
}

// making reisdual graph by overwriting the flow mat,flow[i][j] = CapacityMat[i][j] - flow[i][j];
void create_Augumented_CapacityMat(std::vector < float >& Wmat, std::vector < std::vector <float> >& CapacityMat)
{
	int i = 0;
	for (i = 1; i <= N; i++)
	{
		if (Wmat[i - 1] > 0)
		{
			CapacityMat[0][i] = 0;
			CapacityMat[i][N + 1] = Wmat[i - 1];
		}

		else
		{
			CapacityMat[0][i] = -Wmat[i - 1];
			CapacityMat[i][N + 1] = 0;
		}

	}
}

void create_Residual_graph(std::vector < std::vector <float> >& CapacityMat, std::vector < std::vector <float> >& flow)
{
	int i = 0, j = 0;
	for (i = 0; i < NS; i++)
	{
		for (int j = 0; j < NS; j++)
		{
			flow[i][j] = CapacityMat[i][j] - flow[i][j];
			//cout<< CapacityMat[i][j]-flow[i][j];
		}
		//cout<<endl;
	}
}

void printMatrix(std::vector <std::vector <int> >  & M, int len) {
	int i, j;
	for (i = 0; i < len; i++) {
		for (j = 0; j < len; j++)
			std::cout << "\t" << M[i][j]; 
		std::cout << "\n";
	}
}

void printMatrix(std::vector<int>  & M, int len) {
	int i;
	for (i = 0; i < len; i++) {
		std::cout << "\t" << M[i];

	}std::cout << "\n";
}


void savedata(std::vector< int > & Mat,int l,float del,std::string s)
{
	int w = int(del * 10);
	std::ofstream fobj(s+"-" + std::to_string(l) + "-" + std::to_string(w) + ".csv");
	//fobj << "iteration : " << l << "\ndelta : " << del << std::endl;
	int j;
	Mat.erase(Mat.begin() + 1);
	for (int i = 0; i < Mat.size(); i++)
	{
		std::cout << Mat[i];
	}

		for (int i = 0; i < N; i++)
		{
			if (i%VER == 0 && i)	fobj << std::endl;
			fobj << Mat[i]<<",";
			
		}
	
	
	
	fobj.clear();
}





void savedata(std::vector< float > & Mat, int l, float del, std::string s)
{
	int w = int(del * 10);
	std::ofstream fobj(s + "-" + std::to_string(l) + "-" + std::to_string(w) + ".csv");
	//fobj << "iteration : " << l << "\ndelta : " << del << std::endl;
	
	
		
		for (int i = 0; i < N; i++)
		{

			if (i%VER == 0 && i)	fobj << std::endl;
			fobj << Mat[i] << ",";

		}
	

	fobj.clear();
}



