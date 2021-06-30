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
		if( property == "RANDOMSEED" ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
		}
	}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;


	// Campionamento diretto

	double S_iniziale = 100.;
	double T = 1.;
	double K = 100.;
	double r = 0.1;
	double sigma = 0.25;
	int M = 100000;
	int N = 100;
	int L = int(M/N);

	double C[N] = {0.};
	double C2[N] = {0.};
	double P[N] = {0.};
	double P2[N] = {0.};
	
	for(int i = 0; i < N; i++){
		double sum_call = 0.;
		double sum_put = 0.;

		for(int j = 0; j < L; j++){
			
			double Z = rnd.Gauss(0., 1.);
			double S_finale = S_iniziale * exp((r-pow(sigma,2.)/2.)*T + sigma*Z*sqrt(T));
			sum_call += exp(-r*T) * max(0., S_finale-K);
			sum_put += exp(-r*T) * max(0., K-S_finale);
		}

		C[i] = sum_call / double(L);
		C2[i] = pow(C[i], 2);
		P[i] = sum_put / double(L);
		P2[i] = pow(P[i], 2);
	}

	MediaBlocchi(N, C, C2, "C_diretto", "C_diretto_err");
	MediaBlocchi(N, P, P2, "P_diretto", "P_diretto_err");

	
	// Campionamento discreto

	for(int i = 0; i < N; i++){
		double sum_call = 0.;
		double sum_put = 0.;
		double S_i = 0.;

		for(int j = 0; j < L; j++){

			double Z = rnd.Gauss(0., 1.);
			S_i = S_iniziale * exp((r-pow(sigma,2.)/2.)*(0.01) + sigma*Z*sqrt(0.01));
			
			for(double t = 0.02; t < T+0.01; t = t+0.01){
				Z = rnd.Gauss(0., 1.);
				S_i = S_i * exp((r-pow(sigma,2.)/2.)*(t-(t-0.01)) + sigma*Z*sqrt(t-(t-0.01)));
			}
			sum_call += exp(-r*T) * max(0., S_i-K);
			sum_put += exp(-r*T) * max(0., K-S_i);
		}

		C[i] = sum_call / double(L);
		C2[i] = pow(C[i], 2.);
		P[i] = sum_put / double(L);
		P2[i] = pow(P[i], 2.);
	}

	MediaBlocchi(N, C, C2, "C_discreto", "C_discreto_err");
	MediaBlocchi(N, P, P2, "P_discreto", "P_discreto_err");

	return 0;
}
