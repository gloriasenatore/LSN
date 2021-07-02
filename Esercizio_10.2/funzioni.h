#include "random.h"
#include <vector>
#include <algorithm>
#include <iomanip>
#include <string>

#ifndef __Funzioni__
#define __Funzioni__

using namespace std;

void ToFile(const char* Filename, double * data, int size);
void ToFileApp(const char* Filename, double * data, int size);
double uncertainty(double ave2, double ave, int n);
void MediaBlocchi(int N, double * valori, double * valori2, const char * filename_valori, const char * filename_err);
void MediaBlocchi(int N, double * valori, const char * filename_valori, const char * filename_err);

class Population {

private: 
  Random rnd;
  int N;			//numero di città in ciascun individuo
  int num;			//numero di individui in ogni generazione, è scelto costante, ossia non varia al progredire delle iterazioni
  vector<vector<double> > C;	//vettore che contiene le coppie di coordinate rappresentanti le città
  vector<vector<double> > pop;   //popolazione presente ad una certa generazione

protected:

public:

  // constructors
  Population(int num_i, vector<vector<double> > cities, int rank);
  // destructor
  ~Population();
  // methods
  vector<vector<double> > GetPopulation();			//permette di ottenere la popolazione
  vector<double> GetIndividual(int ind);			//permette di ottenere l'individuo numero ind
  vector<double> PairPermutation(vector<double> G);		//operatore che muta le città a coppie
  bool CheckFunction(vector<double> P);				//funzione di controllo
  void StartingPopulation(vector<double> G);			//imposta la popolazione iniziale partendo da un certo individuo e ottenendone altri attraverso PairPermutation
  double Fitness(vector<double> P);				//permette di ottenere la lunghezza del percorso di un certo individuo
  vector<double> Inversion(vector<double> G);			//operatore che inverte l'ordine di m città contigue
  vector<double> Permutation(vector<double> G);			//operatore che scambia tra loro due gruppi di m città contigue
  void SortingPopulation();					//ordina la popolazione per valori crescenti di lunghezza del percorso
  void BestPath(int rank);					//stampa su file il percorso migliore
  void PrintL(vector<double> P, int rank);			//stampa su file la lunghezza di un percorso
  void PrintAveL(int num_best);					//stampa su file la lunghezza media dei primi num_best individui della popolazione
  int Selection();						//operatore che seleziona gli individui per il crossover
  void Crossover();						//operatore di crossover
  void RandomSearch(int rank);					//metodo per operare una iterazione di random search senza crossover
  void Iteration(int rank);					//metodo per operare una iterazione con crossover

  vector<double> Move(vector<double> G, double beta); 		//metodo che performa l'algoritmo di Metropolis

  void SetIndividual(vector<double> P, int ind);		//setta dall'esterno l'individuo della popolazione al posto ind
};

#endif
