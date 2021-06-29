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

   
	// DISTRIBUZIONE UNIFORME
	
	int N = 100;
	int M = 1000000;
	int L = int(M/N);
	double I[N] = {0.};
	double I2[N] = {0.};

	for(int j = 0; j < N; j++){
		double sum = 0.;

		for(int i = 0; i < L; i++){
			double x = rnd.Rannyu();
			sum += (M_PI / 2.) * cos(M_PI * x / 2.);
		}
	
		I[j] = sum / double(L);
		I2[j] = pow(I[j], 2.);		
	}

	MediaBlocchi(N, I, I2, "2.1_I_uniform", "2.1_I_uniform_err");

	
	// IMPORTANCE SAMPLING
	// Campiono la densità di probabilità con reject

	for(int j = 0; j < N; j++){
		double sum = 0.;

		for(int i = 0; i < L; i++){ 
			int n = 0;

			while(n == 0){
				double x = rnd.Rannyu();
				double r = rnd.Rannyu(); 
				if((3./2.) * r < (3./2.) * (1. - pow(x,2.))){
					sum += (M_PI/2.) * cos(M_PI/2.*x) * (2./3.) / (1 - pow(x,2.));
					n = 1;
				}
			}
		}

		I[j] = sum / double(L);
		I2[j] = pow(I[j], 2.);			
	}

	MediaBlocchi(N, I, I2, "2.1_I_sam", "2.1_I_sam_err");

	return 0;
}
