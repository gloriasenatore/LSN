#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "funzioni.h"
#include "random.h"
#include <math.h>

double a_bohr = 0.0529E-9;

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

double r_sferiche(double x, double y, double z){	
	return sqrt(pow(x,2.) + pow(y,2.) + pow(z,2.));
}

double theta_sferiche(double x, double y, double z){ //theta tra -pi e pi.	
	return acos(z / r_sferiche(x,y,z));
}

double phi_sferiche(double x, double y){ //phi tra 0 e 2pi.	
	return atan2(y,x);
}

double psi_100(double x, double y, double z){
	return pow(a_bohr, -3./2.)/sqrt(M_PI) * exp(-r_sferiche(x,y,z)/a_bohr);
}

double psi_210(double x, double y, double z){
	return pow(a_bohr, -5./2.)/8. * sqrt(2./M_PI) * r_sferiche(x,y,z) * exp(-r_sferiche(x,y,z)/(2.*a_bohr)) * cos(theta_sferiche(x,y,z));
}
