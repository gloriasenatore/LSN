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

	
	//ESERCIZIO 8.2
	//Scrittura della "griglia" di valori di sigma e mu
	bool griglia = true; //set griglia true per valutare l'energia locale dello stato in funzione di sigma e mu
	
	if(griglia == true){
		
		int M = 1000000;
		int N = 100;
		int L = int(M/N);
	
		ofstream griglia;
		griglia.open("8.2.griglia.dat", ios::app);
	
		for(double sigma = 0.58; sigma < 0.64; sigma = sigma+0.01){
			for(double mu = 0.77; mu < 0.83; mu = mu+0.01){

				double traj_iniziale = 0.;
				double H_medio[N] = {0.};
	
				double sum_prog = 0.;
				double su2_prog = 0.;

				for(int i = 0; i < N; ++i){
					double traj = 0.;
					double sum_E = 0.;

					for(int j = 0; j < L; ++j){
			
						double x = rnd.Rannyu(-2.5, 2.5); // campiono con probabilità di transizione T uniforme

						traj = traj_iniziale + x;
		
						double A = min( 1., pow(psi_T(sigma, mu, traj), 2.) / pow(psi_T(sigma, mu, traj_iniziale), 2.) );

						double t = rnd.Rannyu();
						if(t <= A) traj_iniziale = traj;
			
						sum_E += -0.5 * d2_psi_T(sigma, mu, traj_iniziale) / psi_T(sigma, mu, traj_iniziale) +  V(traj_iniziale);
					}

					H_medio[i] = sum_E / double(L);
					sum_prog += H_medio[i];
					su2_prog += H_medio[i]*H_medio[i];
				}

			
				double sum_err = uncertainty(su2_prog/double(N), sum_prog/double(N), N-1);
				griglia << "sigma" << setw(12) << sigma << setw(12) << "mu" << setw(12) << mu << setw(20) << sum_prog/double(N) << setw(20) << sum_err << endl;
			}
		}

		griglia.close();
	}


	//Scelgo sigma=0.62, mu=0.80
	double sigma = 0.62;
	double mu = 0.81;

	ofstream histo;
	histo.open("8.2.histo.dat", ios::app);

	for(double x = -3.; x < 3.; x=x+0.1){
		double psi_2 = pow(psi_T(sigma, mu, x+0.05),2.);
		
		histo << setw(12) << x+0.05 << setw(20) << psi_2 << endl;
	}

	histo.close();


	return 0;
}

