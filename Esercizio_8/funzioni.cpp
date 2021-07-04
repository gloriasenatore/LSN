#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "funzioni.h"
#include "random.h"

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
		return sqrt((ave2 - pow(ave, 2.)) / double(n));
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
			su2_prog[i] += pow(valori[j], 2.);
		}

		sum_prog[i] /= double(i+1); //media cumulativa
		su2_prog[i] /= double(i+1); //media quadra cumulativa

		err[i] = uncertainty(su2_prog[i], sum_prog[i], i);
	}

	ToFile(filename_valori, sum_prog, N);
	ToFile(filename_err, err, N);
}

double psi_T(double sigma, double mu, double x){

	return exp(-pow((x - mu),2.)/(2. * pow(sigma,2.))) + exp(-pow((x + mu),2.)/(2. * pow(sigma,2.)));
}

double d2_psi_T(double sigma, double mu, double x){

	return (1 / pow(sigma,4)) * exp(-pow((x + mu),2.)/(2. * pow(sigma,2.))) * (exp(2.*x*mu/pow(sigma,2.)) * (mu - sigma - x) * (mu + sigma - x) + pow((mu + x),2.) - pow(sigma,2.));
}

double V(double x){

	return pow(x,4.) - 5./2.*pow(x,2.);
}
