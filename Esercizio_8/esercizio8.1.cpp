#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "funzioni.h"
#include "random.h"

using namespace std;

int main(int argc, char *argv[]){

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

	//ESERCIZIO 8.1
	int M = 1000000;
	int N = 100;
	int L = int(M/N);

	double traj_iniziale = 0.;
	double H_medio[N] = {0.};

	double sigma = 0.60;
	double mu = 0.80;

	for(int i = 0; i < N; ++i){
		double traj = 0.;
		double sum_E = 0.;
		double sum_A = 0.;

		for(int j = 0; j < L; ++j){
			
			double x = rnd.Rannyu(-2.5, 2.5); // campiono con probabilità di transizione T uniforme

			traj = traj_iniziale + x;
		
			double A = min( 1., pow(psi_T(sigma, mu, traj), 2.) / pow(psi_T(sigma, mu, traj_iniziale), 2.) );
			sum_A += A;

			double t = rnd.Rannyu();
			if(t <= A) traj_iniziale = traj;
			
			sum_E += -0.5 * d2_psi_T(sigma, mu, traj_iniziale) / psi_T(sigma, mu, traj_iniziale) +  V(traj_iniziale);
		}

		H_medio[i] = sum_E / double(L);
		cout << "accettanza media = " << sum_A/(double(L)) << endl;
	}

	MediaBlocchi(N, H_medio, "8.1.H_medio", "8.1.H_err");

	return 0;
}
