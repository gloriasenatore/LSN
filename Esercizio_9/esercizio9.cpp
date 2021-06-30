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

	//32 CITIES ON A CIRCUMFERENCE ---> set square false
	//32 CITIES ON A SQUARE        ---> set square true  
	bool square = false;
	if(square == true) cout << "Problema del commesso viaggiatore. 32 città sulla superficie di un quadrato" << endl;
	if(square == false) cout << "Problema del commesso viaggiatore. 32 città sulla circonferenza" << endl;

	//Creation of the first gene
	vector<vector<double>> C;	//Vettore che contiene le città, ossia le coppie di coordinate
	ofstream cities;
	cities.open("cities.dat");

	double x = 0.;
	double y = 0.;
	for(int i = 0; i < 32; ++i){
		if(square == false){
			x = rnd.Rannyu(-1., 1.);
			y = 0.;
			double f = rnd.Rannyu();
			if(f < 0.5) y = - sqrt(1. - x*x);	//Vincoli della circonferenza
			else y = sqrt(1. - x*x);
		}

		else{
			x = rnd.Rannyu(-1., 1.);       //Vincoli della superficie quadrata
			y = rnd.Rannyu(-1., 1.);
		}

		vector<double> v = {x, y};
				
		C.push_back(v);
	}

	vector<double> G(32);
	for(int i = 0; i < 32; ++i) {
		cities << setw(20) << C[i][0] << setw(20) << C[i][1] << endl;
		G[i] = (double)i;
	}
	cities.close();

	Population Genetic(100, C);	//creo la classe "Genetic" inizializzandola con 100 individui, per città sulla circonferenza, e 200, per città nel quadrato, e con il vettore delle città

	Genetic.StartingPopulation(G);

	Genetic.SortingPopulation();

	vector<vector<double>> P;
	P = Genetic.GetPopulation();
	Genetic.PrintL(P[0]);

	//RANDOM SEARCH	---> prima di includere il crossover, verifico che il GA performi un buon random search con i soli operatori di mutazione. Set randomsearch true per fare ciò
	
	bool randomsearch = false;

	if(randomsearch == true){

		cout << "Random search con i soli operatori di mutazione" << endl;
		for(int i = 0; i < 60000; ++i){
			if(i%1000 == 0) cout << "iteration number " << i << endl;
 			Genetic.RandomSearch();
		}	

		Genetic.BestPath();
	}

	//CROSSOVER

	else{
	
		cout << "Algoritmo genetico con operatori di mutazione e di crossover" << endl;
		for(int i = 0; i < 20000; ++i){
			if(i%1000 == 0) cout << "iteration number " << i << endl;
			Genetic.Iteration();
		}
		Genetic.BestPath();
	}

	return 0;
}
