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

	// VALOR MEDIO

	int N = 100;
	int M = 100000;
	int L = int(M/N);
	
	double ave[N] = {0.};
	double av2[N] = {0.};

	for(int i = 0; i < N; i++){
		double sum = 0.;

		for(int j = 0; j < L; j++){
			sum += rnd.Rannyu();
		}

		ave[i] = sum / double(L); //valor medio di ogni blocco
		av2[i] = (ave[i]) * (ave[i]); //valor medio quadro di ogni blocco
	}

	MediaBlocchi(N, ave, av2, "Valor_medio_1.1", "Valor_medio_err_1.1");

	// VARIANZA

	double avev[N] = {0.};
	double av2v[N] = {0.};

	for(int i = 0; i < N; i++){
		double sum = 0.;

		for(int j = 0; j < L; j++){
			sum += pow((rnd.Rannyu() - 0.5), 2.);
		}

		avev[i] = sum / double(L); //valor medio di ogni blocco
		av2v[i] = (avev[i]) * (avev[i]); //valor medio quadro di ogni blocco
	}

	MediaBlocchi(N, avev, av2v, "Valor_medio_var_1.1", "Valor_medio_var_err_1.1");
	

	// TEST DEL CHI QUADRO
	
	int m = 100; // numero di sotto-intervalli
	int num = 10000; // numero di numeri casuali
	int n[m] = {0}; // array i cui elementi sono il numero di numeri casuali contenuti in ogni sottointervallo
	double chi[100] = {0.};

	for(int i = 0; i < 100; i++){
		double chi2 = 0.;
		for(int i = 0; i < 100; i++){
			n[i] = 0;
		}
		
		for(int l = 0; l < num; l++){		
			int random_number = int(rnd.Rannyu(0., 100.));
			
			for(int k = 0; k < m; k++){
				if (random_number == k){
					n[k] += 1;
				}
			} 
		}

		for(int l = 0; l < m; l++){
			chi2 += pow((double(n[l]) - (double(num) / double(m))), 2.) / (double(num) / double(m));
		}
		chi[i] = chi2;
	}

	ToFile("Chi2_1.1", chi, 100);

	return 0;
}
