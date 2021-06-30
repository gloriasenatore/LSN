#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "funzioni.h"

using namespace std;

void ToFile(const char* Filename, double * data, int size){

	ofstream fout(Filename);
	if (fout.is_open()){
		for (int i = 0; i < size; i++) fout << data[i] << endl;
	}
	else cerr << "PROBLEM: Unable to open " << Filename << endl;

	fout.close();
}

void ToFileApp(const char* Filename, double * data, int size){

	ofstream fout(Filename, ios::app);
	if (fout.is_open()){
		for (int i = 0; i < size; i++) fout << data[i] << endl;
	}
	else cerr << "PROBLEM: Unable to open " << Filename << endl;

	fout.close();
}

double uncertainty(double ave2, double ave, int n){

	if(n == 0){
		return 0;
	}

	else{
		return sqrt((ave2 - pow(ave, 2)) / double(n));
	}
}

void MediaBlocchi(int N, double * valori, double * valori2, const char * filename_valori, const char * filename_err){

	double sum_prog[N] = {0.};
	double su2_prog[N] = {0.};
	double err[N] = {0.};

	for(int i = 0; i < N; i++){
		for(int j = 0; j < i+1; j++){
			sum_prog[i] += valori[j];
			su2_prog[i] += valori2[j];
		}

		sum_prog[i] /= double(i+1); //media cumulativa
		su2_prog[i] /= double(i+1); //media quadra cumulativa

		err[i] = uncertainty(su2_prog[i], sum_prog[i], i);
	}

	ToFile(filename_valori, sum_prog, N);
	ToFile(filename_err, err, N);
}

void MediaBlocchi(int N, double * valori, const char * filename_valori, const char * filename_err){

	double sum_prog[N] = {0.};
	double su2_prog[N] = {0.};
	double err[N] = {0.};

	for(int i = 0; i < N; i++){
		for(int j = 0; j < i+1; j++){
			sum_prog[i] += valori[j];
			su2_prog[i] += pow(valori[j], 2);
		}

		sum_prog[i] /= double(i+1); //media cumulativa
		su2_prog[i] /= double(i+1); //media quadra cumulativa

		err[i] = uncertainty(su2_prog[i], sum_prog[i], i);
	}

	ToFile(filename_valori, sum_prog, N);
	ToFile(filename_err, err, N);
}

Population :: Population(int num_i, vector<vector<double>> cities){

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

	N = cities.size();
	num = num_i;
	C = cities;
}

Population :: ~Population(){}

vector<vector<double>> Population :: GetPopulation(){

	return pop;
}

vector<double> Population :: PairPermutation(vector<double> G){

	vector<double> P(N);
	for(int i = 0; i < N; ++i) P[i] = G[i];

	int f = (int)rnd.Rannyu(1., double(N));
	int g = (int)rnd.Rannyu(1., double(N));
	while(g == f){
		f = (int)rnd.Rannyu(1., double(N));	//Non devo selezionare la stessa città, né la prima
	}

	P[g] = G[f];
	P[f] = G[g];

	ofstream problems;
	problems.open("problems.dat",ios::app);
	bool check = CheckFunction(P);
	if(check == false) problems << "Problemi con PairPermutation" << endl;
	problems.close();
		
	return P;

}

vector<double> Population :: Permutation(vector<double> G){

	G.erase(G.begin());	//Cancello la prima città, che non deve permutare
	
	vector<double> P(N-1);
	for(int i = 0; i < N-1; ++i) P[i] = G[i];

	int m = (int)rnd.Rannyu(1., double(N-1)/2.);			//numero città contigue da permutare con altrettante città
	int c1 = (int)rnd.Rannyu(0., double(N-1)/2. - double(m));	//prima città del primo gruppo, scelta nella prima metà dell'individuo
	int c2 = (int)rnd.Rannyu(double(c1+m), double(N-1-m));		//prima città del secondo gruppo, scelta nella seconda metà dell'individuo
	
	if(c1 == c2) c2 = c2+1;	//evito che i due gruppi da permutare si sovrappongano

	for(int i = 0; i < m; ++i){
		P[c1+i] = G[c2+i];
		P[c2+i] = G[c1+i];
	}

	G.insert(G.begin(), 0);		//riaggiungo la prima città a G e a P
	P.insert(P.begin(), 0);

	ofstream problems;
	problems.open("problems.dat",ios::app);
	bool check = CheckFunction(P);
	if(check == false) problems << "Problemi con Permutation" << endl;
	problems.close();

	return P;
}

vector<double> Population :: Inversion(vector<double> G){

	vector<double> P(N);
	for(int i = 0; i < N; ++i) P[i] = G[i];	

	int m = (int)rnd.Rannyu(2., double(N));		//numero città contigue il cui ordine è da invertire
	int c1 = (int)rnd.Rannyu(1., double(N-m));	//prima città da cui partire per invertire

	for(int i = 0; i < m; ++i){
		P[c1+i] = G[c1+m-1-i];
	}

	ofstream problems;
	problems.open("problems.dat",ios::app);
	bool check = CheckFunction(P);
	if(check == false) problems << "Problemi con Inversion" << endl;
	problems.close();

	return P;
}

bool Population :: CheckFunction(vector<double> P){

	for(int i = 0; i < N; ++i){
		for(int j = i+1; j < N; ++j){
			if(P[i] == P[j]){
				cout << "Problema: la città " << P[i] << " viene ripetuta ai posti " << i << " e " << j << endl;
				return false;
			}
		} 
	}

	if(P[0] != 0.){
		cout << "Attenzione: problemi con la prima città. Ricontrollare il codice" << endl;
		return false;  
	}

	return true;
}

int Population :: Selection(){

	double r = rnd.Rannyu();
	int j = int(num * pow(r, 10.));
	return j;
}

void Population :: StartingPopulation(vector<double> G){

	pop.push_back(G);

	vector<double> T(N);
	for(int i = 0; i < N; ++i) T[i] = G[i]; 

	for(int n = 1; n < num; ++n){
		vector<double> P(N);
		P = PairPermutation(T);
		pop.push_back(P);
		for(int i = 0; i < N; ++i) T[i] = P[i];	//ogni nuovo individuo è ottenuto permutando una coppia di città nell'individuo precedente
		
	}	
}

double Population :: Fitness(vector<double> P){

	double L = 0.;
	
	for(int i = 0; i < N-1; ++i){
		int l1 = P[i];	//l'ordine delle città lo trovo nel "genoma" degli individui
		int l2 = P[i+1];
		L += pow(C[l1][0] - C[l2][0], 2.) + pow(C[l1][1] - C[l2][1], 2.);	//norma L2
	}

	int l1 = P[31];
	int l2 = P[0];
	L += pow(C[l1][0] - C[l2][0], 2.) + pow(C[l1][1] - C[l2][1], 2.);
	
	return L;
}

void Population :: SortingPopulation(){

	vector<double> L(num);
	for(int n = 0; n < num; ++n){
		L[n] = Fitness(pop[n]);
		pop[n][N] = L[n];	//aggiungo a ciascun individuo l'informazione della sua lunghezza
	}

	sort(pop.begin(), pop.end(), [&](const vector<double>& a, const vector<double>& b) {return a[N] < b[N];});	//Riordino gli elementi in modo crescente di L
}

void Population :: BestPath(){

	ofstream path;
	path.open("bestpath.dat");

	for(int i = 0; i < N; ++i){
		int l = pop[0][i];
		if(i == 0) path << setw(20) << C[l][0] << setw(20) << C[l][1] << setw(20) << l << setw(20) << Fitness(pop[0]) << endl;	
		else path << setw(20) << C[l][0] << setw(20) << C[l][1] << setw(20) << l << endl;
	}
	path.close();
}

void Population :: PrintL(vector<double> P){

	ofstream L;
	L.open("bestL.dat",ios::app);

	L << Fitness(P) << endl;
	
	L.close();
}

void Population :: PrintAveL(int num_best){

	ofstream L;
	L.open("bestAveL.dat",ios::app);

	double fit = 0.;
	for(int i = 0; i < num_best; ++i){
		fit += Fitness(pop[i]);
	}
	
	L << fit / double(num_best) << endl;
	L.close();
}

void Population :: Crossover(){

	int mother = Selection();
	int father = Selection();
	while(father == mother) father = Selection();	//padre e madre devono essere due individui diversi
	int f = rnd.Rannyu(1, N-3);	//posizione in cui vengono "tagliati" i cromosomi dei due individui scelti

	vector<double> son1(f+1);
	vector<double> son2(f+1);
	for(int i = 0; i <= f; ++i){
		son1[i] = pop[mother][i];	//la prima parte del cromosoma dei due figli è uguale a quella di madre e padre rispettivamente
		son2[i] = pop[father][i];
	}

	for(int i = f+1; i < N; ++i){
		int p = 1;
		while(find(son1.begin(), son1.end(), pop[father][p]) != son1.end()) p++;	//la seconda parte contiene gli elementi dell'altro genitore, non presenti nella prima parte, nell'ordine in cui appaiono nel genitore
		son1.push_back(pop[father][p]);

		int q = 1;
		while(find(son2.begin(), son2.end(), pop[mother][q]) != son2.end()) q++; 
		son2.push_back(pop[mother][q]);		
	}
	
	bool check = CheckFunction(son1);
	ofstream problems;
	problems.open("problems.dat",ios::app);
	if(check == false) problems << "Problemi con il crossover per figlio1" << endl;
	check = CheckFunction(son2);
	if(check == false) problems << "Problemi con il crossover per figlio2" << endl;

	problems.close();

	for(int i = 0; i < N; ++i) pop[num-2][i] = son1[i];	//sostituisco i figli agli ultimi due individui della poplazione ordinata
	
	for(int i = 0; i < N; ++i) pop[num-1][i] = son2[i];
}

void Population :: RandomSearch(){

	double f = rnd.Rannyu();
	int c = (int)rnd.Rannyu(0., num);
	vector<double> P(N);

	if(f < 0.4){ 
		P = PairPermutation(pop[c]);
		for(int i = 0; i < N; ++i) pop[num-1][i] = P[i];
	}
	if(0.4 < f && f < 0.55){ 
		P = Permutation(pop[c]);
		for(int i = 0; i < N; ++i) pop[num-1][i] = P[i];
	}
	if(0.55 < f && f < 0.8){
		P = Inversion(pop[c]);
		for(int i = 0; i < N; ++i) pop[num-1][i] = P[i];
	}

	SortingPopulation();
	PrintL(pop[0]);
}

void Population :: Iteration(){

	double pc = rnd.Rannyu();	//con probabilità del 60% avviene il crossover
	if(pc < 0.6){
		Crossover();
			
		double pm = rnd.Rannyu();

		if(pm < 0.1){	//mutazioni sui figli
			vector<double> P(N);
			vector<double> Q(N);
			P = PairPermutation(pop[num-2]);
			Q = PairPermutation(pop[num-1]);
			for(int t = 0; t < N; ++t){
				pop[num-2][t] = P[t];	//non avendo già ordinato la popolazione in Crossover(), sono sicura che i figli saranno al penultimo e all'ultimo posto della popolazione
				pop[num-1][t] = Q[t];
			}
		}
		if(0.1 < pm && pm < 0.15){ 
			vector<double> P(N);
			vector<double> Q(N);
			P = Permutation(pop[num-2]);
			Q = Permutation(pop[num-1]);
			for(int t = 0; t < N; ++t){
				pop[num-2][t] = P[t];
				pop[num-1][t] = Q[t];
			}
		}
		if(0.2 < pm && pm < 0.3){
			vector<double> P(N);
			vector<double> Q(N);
			P = Inversion(pop[num-2]);
			Q = Inversion(pop[num-1]);
			for(int t = 0; t < N; ++t){
				pop[num-2][t] = P[t];
				pop[num-1][t] = Q[t];
			}
		}
		
		
		SortingPopulation();
	}
	PrintL(pop[0]);
	PrintAveL(num/2.);	//stampo la lunghezza media della prima metà della popolazione ordinata
}
