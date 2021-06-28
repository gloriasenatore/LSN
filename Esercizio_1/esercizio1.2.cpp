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


	int N[4] = {1, 2, 10, 100};
		
	for(int i = 0; i < 4; i++){
		double std[10000] = {0.};
		double exp[10000] = {0.};
		double lor[10000] = {0.};

		for(int j = 0; j < 10000; j++){
			double sum_std = 0.;
			double sum_exp = 0.;
			double sum_lor = 0.;
			
			for(int k = 0; k < N[i]; k++){
				sum_std += rnd.Rannyu();
				sum_exp += rnd.Exponential(1.);
				sum_lor += rnd.Lorentzian(1., 0.);
			}

			std[j] = sum_std / double(N[i]);
			exp[j] = sum_exp / double(N[i]);
			lor[j] = sum_lor / double(N[i]);
		}

		ToFileApp("Risultati_std_1.2", std, 10000);
		ToFileApp("Risultati_exp_1.2", exp, 10000);
		ToFileApp("Risultati_lor_1.2", lor, 10000);
	}


	return 0;
}
