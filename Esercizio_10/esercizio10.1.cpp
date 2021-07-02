#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "funzioni.h"
#include <cmath>
#include <vector>
#include <iomanip>

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

	
	//SIMULATED ANNEALING
	vector<double> x(32);
	vector<double> y(32);
	fstream cities;
	cities.open("cities.dat");	//utilizzo lo stesso vettore di città dell'esercizio 9 per poter confrontare i risultati
	for(int i = 0; i < 32; ++i){		
		cities >> x[i] >> y[i];
	}

	cities.close();

	vector<double> G(32);
	vector<vector<double>> C;
	for(int i = 0; i < 32; ++i){
		vector<double> v = {x[i], y[i]};
		C.push_back(v);
		G[i] = (double)i;
	}

	Population SA(1, C);	//creo la classe "SA" inizializzandola con un individuo e con il vettore delle città
	SA.StartingPopulation(G);
	SA.PrintL(G);

	cout << "Simulated annealing per risolvere l'algoritmo del commesso viaggiatore" << endl;

	vector<double> beta = {10., 100., 1000., 10000};
	for(int i = 1; i <= 40000; ++i){
		if(i%5000 == 0) cout << "iteration number " << i << endl;
		if(i <= 10000){
			G = SA.Move(G, beta[0]);
			SA.PrintL(G);
		}
		if(10000 < i && i  <= 20000){
			G = SA.Move(G, beta[1]);
			SA.PrintL(G);
		}
		if(20000 < i && i <= 30000){
			G = SA.Move(G, beta[2]);
			SA.PrintL(G);
		}
		if(30000 < i && i <= 40000){
			G = SA.Move(G, beta[3]);
			SA.PrintL(G);
		}
	 }
	SA.BestPath();

	return 0;
}
