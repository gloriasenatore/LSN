#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "funzioni.h"
#include <math.h>

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
		if( property == "RANDOMSEED" ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
		}
	}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	double L = 0.9;
	double d = 1.;

	int N = 100;
	int M = 10000000;
	int F = int(M/N); //numero tiri in ogni blocco
	double pi_estimation[N] = {0.};
	double x = 0.;
	double theta = 0.;


	// Metodo ricorsivo

	for(int i = 0; i < N; i++){
		int N_hit = 0;
		//double p = 2.;
		//double p = 0.1;
		double p = 3.;

		for(int j = 0; j < F; j++){
			x = rnd.Rannyu(0., d/2.);
			theta = rnd.Rannyu(0., p/2.);

			if(x <= (L/2.) * sin(theta)){
				N_hit++;
			}
			if(N_hit != 0){
				p = ((2. * L * double(j)) / (double(N_hit) * d));
			}		
		}

		pi_estimation[i] = ((2. * L * double(F)) / (double(N_hit) * d));
	}

	//MediaBlocchi(N, pi_estimation, "Somme_progressive_p", "Err_p");
	//MediaBlocchi(N, pi_estimation, "Somme_progressive_p_0.1", "Err_p_0.1");
	MediaBlocchi(N, pi_estimation, "Somme_progressive_p_3", "Err_p_3");

	// Metodo reject

	for(int i = 0; i < N; i++){
		int N_hit = 0;
		double theta = 0.;

		for(int j = 0; j < F; j++){
			int n = 0;
			while(n == 0){
				double x = rnd.Rannyu(0., 1.);
				double y = rnd.Rannyu(0., 1.);

				if(pow(x,2)+pow(y,2.) < 1){
					theta = acos(x / sqrt(pow(x,2.)+pow(y,2.)));
					n++;
				}
			}
			double r = rnd.Rannyu(0., d/2.);
			if(r <= (L/2.) * sin(theta)){
				N_hit++;
			}		
		}

		pi_estimation[i] = ((2. * L * double(F)) / (double(N_hit) * d));
	}

	MediaBlocchi(N, pi_estimation, "Somme_progressive_p_rej", "Err_p_rej");

	return 0;
}
