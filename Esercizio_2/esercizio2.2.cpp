#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "funzioni.h"
#include <cmath>

using namespace std;
 
int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	fstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
		if( property == "SEED2" ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
		}
	}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;


	// DISCRETE LATTICE

	double squared_displacement[10000][101] = { {0., 0.} };

	for(int n = 0; n < 10000; n++){

		int traj[101][3] = { {0, 0} };

		for(int i = 0; i < 101; i++){

			if(i != 0){
				int u = int(rnd.Rannyu(0., 3.));
				double p = rnd.Rannyu();
				for(int j = 0; j < 3; j++){
					traj[i][j] = traj[i-1][j];
				}

				if(p < 0.5){
					traj[i][u] = traj[i-1][u] - 1;
				}
				else{
					traj[i][u] = traj[i-1][u] + 1;
				}
			}
			
			squared_displacement[n][i] = pow(traj[i][0],2.) + pow(traj[i][1],2.) + pow(traj[i][2],2.);
		}
	}

	double mean_squared_displacement[101] = {0.};
	double sqrt_msd[101] = {0.};
	double mean_squared_displacement2[101] = {0.};
	double sqrt_msd2[101] = {0.};
	double err[101] = {0.};
	
	for(int i = 0; i < 101; i++){
		for(int n = 0; n < 10000; n++){
			mean_squared_displacement[i] += squared_displacement[n][i];
			mean_squared_displacement2[i] += pow(squared_displacement[n][i], 2.);
		}
		sqrt_msd[i] = sqrt(double(mean_squared_displacement[i]) / double(10000));
		sqrt_msd2[i] = sqrt(double(mean_squared_displacement2[i]) / double(10000));
		err[i] = uncertainty(sqrt_msd2[i], sqrt_msd[i], 10000-1);
	}
	
	ToFile("2.2_sqrtmean_discrete", sqrt_msd, 101);
	ToFile("2.2_err_discrete", err, 101);
	
	// CONTINUE LATTICE

	for(int n = 0; n < 10000; n++){

		double traj[101][3] = { {0, 0} };

		for(int i = 0; i < 101; i++){

			if(i != 0){
				double phi = rnd.Rannyu(0., 2.*M_PI);
				double r = rnd.Rannyu();
				double theta = acos(1. - 2.*r);

				for(int j = 0; j < 3; j++){
					traj[i][j] = traj[i-1][j];
				}

				traj[i][0] = traj[i-1][0] + sin(theta)*cos(phi);
				traj[i][1] = traj[i-1][1] + sin(theta)*sin(phi);
				traj[i][2] = traj[i-1][2] + cos(theta);
			}

			squared_displacement[n][i] = pow(traj[i][0],2) + pow(traj[i][1],2) + pow(traj[i][2],2);
		}
	}
	
	for(int i = 0; i < 101; i++){
		mean_squared_displacement[i] = 0;
		mean_squared_displacement2[i] = 0;
		for(int n = 0; n < 10000; n++){
			mean_squared_displacement[i] += squared_displacement[n][i];
			mean_squared_displacement2[i] += pow(squared_displacement[n][i], 2.);
		}
		sqrt_msd[i] = sqrt(mean_squared_displacement[i] / double(10000));
		sqrt_msd2[i] = sqrt(double(mean_squared_displacement2[i]) / double(10000));
		err[i] = uncertainty(sqrt_msd2[i], sqrt_msd[i], 10000-1);
	
	}
	
	ToFile("2.2_sqrtmean_continue", sqrt_msd, 101);
	ToFile("2.2_err_continue", err, 101);

	return 0;
}
