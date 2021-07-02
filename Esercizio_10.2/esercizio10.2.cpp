#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "funzioni.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "mpi.h"

using namespace std;
 
int main (int argc, char *argv[]){

	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double tstart = MPI_Wtime(); //Valuto le performance dell'algoritmo
	if(size != 4) { cout << "Usare quattro core" << endl; return 1; }

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
			input >> property;	//Setto un seme diverso per ogni core
			if(rank == 0 && property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
			if(rank == 1 && property == "SEED2" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
			if(rank == 2 && property == "SEED3" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
			if(rank == 3 && property == "SEED4" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
	
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;


	//32 CITIES ON A SQUARE 
	if(rank==0) cout << "Problema del commesso viaggiatore con programmazione parallela. 32 città sulla superficie di un quadrato. Programmazione parallela" << endl;

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

	Population Genetic(200, C, rank);	//creo la classe "Genetic" inizializzandola con 200 individui e con il vettore delle città

	Genetic.StartingPopulation(G);
	Genetic.SortingPopulation();

	vector<vector<double> > P;
	P = Genetic.GetPopulation();
	Genetic.PrintL(P[0], rank);


	//CROSSOVER

	MPI_Status stat1, stat2, stat3, stat4;
	MPI_Request req1, req2, req3, req4;
	int itag = 1; int itag2 = 2; int itag3 = 3; int itag4 = 4;

	int Nmigr = 200; //Numero di generazioni tra una migrazione e l'altra

	
	
	if(rank==0) cout << "Algoritmo genetico con operatori di mutazione e di crossover. Migrazione degli individui migliori tra i contentinenti" << endl;
	for(int i = 0; i < 20000; ++i){
		if(i%1000 == 0 && rank == 0) cout << "iteration number " << i << endl;
		Genetic.Iteration(rank);

		if(i%Nmigr == 0){
			MPI_Barrier(MPI_COMM_WORLD); // Aspetto che tutti i processi siano arrivati a questo punto, cioè abbiano compiuto n.intero*Nmigr iterazioni

			vector<int> numbers = {0, 1, 2, 3};
			vector <int> Crecv(4);
			if(rank == 0) {
			 Crecv[0] = int(rnd.Rannyu(1., 4.)); //Crecv[i] è il rank che riceve l'individuo dal rank i-esimo. Deve essere diverso dal rank da cui 													riceve e dagli altri Crecv[i].
			 remove(numbers.begin(), numbers.end(), Crecv[0]);

			 Crecv[1] = numbers[int(rnd.Rannyu(0., 3.))]; 
			 while(Crecv[1] == 1) Crecv[1] = numbers[int(rnd.Rannyu(0., 3.))];
			 remove(numbers.begin(), numbers.end(), Crecv[1]);
	
			 Crecv[2] = numbers[int(rnd.Rannyu(0., 2.))]; 
 			 while(Crecv[2] == 2) Crecv[2] = numbers[int(rnd.Rannyu(0., 2.))];
			 remove(numbers.begin(), numbers.end(), Crecv[2]);

			 Crecv[3] = numbers[0];
			 while(Crecv[3] == 3){
			  int temp = Crecv[3];
			  Crecv[3] = Crecv[2];
			  Crecv[2] = temp;
			 }
			}

			MPI_Bcast(&Crecv[0], 4, MPI_INTEGER, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD); // Aspetto che tutti i processi abbiano il rank a cui cedere l'individuo

			vector <double> A(32); vector <double> B(32); vector <double> C(32); vector <double> D(32); 
			if(rank == 0) {A = Genetic.GetIndividual(0); MPI_Isend(&A[0], 32, MPI_DOUBLE_PRECISION, Crecv[0], itag, MPI_COMM_WORLD, &req1);}
			if(rank == Crecv[0]) MPI_Recv(&A[0], 32, MPI_DOUBLE_PRECISION, 0, itag, MPI_COMM_WORLD, &stat1);

			if(rank == 1) {B = Genetic.GetIndividual(0); MPI_Isend(&B[0], 32, MPI_DOUBLE_PRECISION, Crecv[1], itag2, MPI_COMM_WORLD, &req2);}
			if(rank == Crecv[1]) MPI_Recv(&B[0], 32, MPI_DOUBLE_PRECISION, 1, itag2, MPI_COMM_WORLD, &stat2);

			if(rank == 2) {C = Genetic.GetIndividual(0); MPI_Isend(&C[0], 32, MPI_DOUBLE_PRECISION, Crecv[2], itag3, MPI_COMM_WORLD, &req3);}
			if(rank == Crecv[2]) MPI_Recv(&C[0], 32, MPI_DOUBLE_PRECISION, 2, itag3, MPI_COMM_WORLD, &stat3);

			if(rank == 3) {D = Genetic.GetIndividual(0); MPI_Isend(&D[0], 32, MPI_DOUBLE_PRECISION, Crecv[3], itag4, MPI_COMM_WORLD, &req4);}
			if(rank == Crecv[3]) MPI_Recv(&D[0], 32, MPI_DOUBLE_PRECISION, 3, itag4, MPI_COMM_WORLD, &stat4); 

			MPI_Barrier(MPI_COMM_WORLD); // Aspetto che tutti i processi si siano scambiati gli individui migliori

			if(rank == Crecv[0]) Genetic.SetIndividual(A, 0); // Ogni continente sostituisce il proprio migliore individuo con quello di un altro continente
			if(rank == Crecv[1]) Genetic.SetIndividual(B, 0);
			if(rank == Crecv[2]) Genetic.SetIndividual(C, 0);
			if(rank == Crecv[3]) Genetic.SetIndividual(D, 0);
				
	
		}
	}

	Genetic.BestPath(rank);
	
	double tend = MPI_Wtime();
	cout << "Tempo impiegato: " << tend - tstart <<  ". Sono il rank numero " << rank << endl;

	MPI_Finalize();

	return 0;
}
