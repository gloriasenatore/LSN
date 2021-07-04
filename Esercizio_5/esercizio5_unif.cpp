#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
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

	int M = 1000000;
	int N = 100;
	int L = int(M/N);

	double a_bohr = 0.0529E-9;
	double traj_iniziale[3] = {a_bohr*3./2., 0., 0.};
	//double traj_iniziale[3] = {a_bohr*150, 0., 0.}; // Scommentare e commentare la parte di equilibrazione per osservare cosa succede con una condizione iniziale per cui l'elettrone è molto distante dal nucleo
	double r_medio[N] = {0.};
	ofstream fout("unif_punti_100", ios::app);
	//ofstream fout("unif_punti_100_lontano", ios::app); // Scommentare per salvare su file i punti partendo dalla condizione iniziale di elettrone molto lontano

	// T UNIFORME
	// psi_100
	cout << "Valor medio di psi_100 con T uniforme " << endl;

	//Equilibrazione del campionamento per 200 passi
	double traj[3] = {0., 0., 0.};
	for(int j = 0; j < 200; ++j){
			
		double phi = rnd.Rannyu(0., 2.*M_PI);
		double f = rnd.Rannyu();
		double theta = acos(1. - 2.*f);
		double r = rnd.Rannyu(0., 2.5*a_bohr);

		traj[0] = traj_iniziale[0] + r*sin(theta)*cos(phi);
		traj[1] = traj_iniziale[1] + r*sin(theta)*sin(phi);
		traj[2] = traj_iniziale[2] + r*cos(theta);
		
		double A = min( 1., pow(psi_100(traj[0], traj[1], traj[2]), 2.) / pow(psi_100(traj_iniziale[0], traj_iniziale[1], traj_iniziale[2]), 2.) );

		double t = rnd.Rannyu();
		if(t <= A){
			for(int n = 0; n < 3; ++n){
				traj_iniziale[n] = traj[n];
			}
		}
	}


	for(int i = 0; i < N; ++i){
		double sum_r = 0.;
		double sum_A = 0.;

		for(int j = 0; j < L; ++j){
			for(int n = 0; n < 3; ++n) fout << traj_iniziale[n]/a_bohr << endl;
			
			double phi = rnd.Rannyu(0., 2.*M_PI);
			double f = rnd.Rannyu();
			double theta = acos(1. - 2.*f);
			double r = rnd.Rannyu(0., 2.5*a_bohr);

			traj[0] = traj_iniziale[0] + r*sin(theta)*cos(phi);
			traj[1] = traj_iniziale[1] + r*sin(theta)*sin(phi);
			traj[2] = traj_iniziale[2] + r*cos(theta);
		
			double A = min( 1., pow(psi_100(traj[0], traj[1], traj[2]), 2.) / pow(psi_100(traj_iniziale[0], traj_iniziale[1], traj_iniziale[2]), 2.) );
			sum_A += A;

			double t = rnd.Rannyu();
			if(t <= A){
				for(int n = 0; n < 3; ++n){
					traj_iniziale[n] = traj[n];
				}
			}
			
			sum_r += r_sferiche(traj_iniziale[0], traj_iniziale[1], traj_iniziale[2])/a_bohr;
		}

		r_medio[i] = sum_r / double(L);
		cout << "accettanza: " << sum_A/(double(L)) << endl;
	}

	fout.close();
	MediaBlocchi(N, r_medio, "unif_r100", "unif_r100_err");
	//MediaBlocchi(N, r_medio, "unif_r100_lontano", "unif_r100_err_lontano");

	// psi_210
	traj_iniziale[0] = 0;
	traj_iniziale[1] = 0;
	traj_iniziale[2] = a_bohr*5.;
	ofstream fout2("unif_punti_210", ios::app);
	cout << "Valor medio di psi_210 con T uniforme " << endl;

	//Equilibrazione del campionamento per 200 passi
	for(int j = 0; j < 200; ++j){
			
		double phi = rnd.Rannyu(0., 2.*M_PI);
		double f = rnd.Rannyu();
		double theta = acos(1. - 2.*f);
		double r = rnd.Rannyu(0., 5.8*a_bohr);

		traj[0] = traj_iniziale[0] + r*sin(theta)*cos(phi);
		traj[1] = traj_iniziale[1] + r*sin(theta)*sin(phi);
		traj[2] = traj_iniziale[2] + r*cos(theta);
		
		double A = min( 1., pow(psi_210(traj[0], traj[1], traj[2]), 2.) / pow(psi_210(traj_iniziale[0], traj_iniziale[1], traj_iniziale[2]), 2.) );

		double t = rnd.Rannyu();
		if(t <= A){
			for(int n = 0; n < 3; ++n){
				traj_iniziale[n] = traj[n];
			}
		}
	}

	for(int i = 0; i < N; ++i){
		double sum_r = 0.;
		double sum_A = 0.;
		r_medio[i] = 0.;

		for(int j = 0; j < L; ++j){
			for(int n = 0; n < 3; ++n) fout2 << traj_iniziale[n]/a_bohr << endl;
			
			double phi = rnd.Rannyu(0., 2.*M_PI);
			double f = rnd.Rannyu();
			double theta = acos(1. - 2.*f);
			double r = rnd.Rannyu(0., 5.8*a_bohr);

			traj[0] = traj_iniziale[0] + r*sin(theta)*cos(phi);
			traj[1] = traj_iniziale[1] + r*sin(theta)*sin(phi);
			traj[2] = traj_iniziale[2] + r*cos(theta);
		
			double A = min( 1., pow(psi_210(traj[0], traj[1], traj[2]), 2.) / pow(psi_210(traj_iniziale[0], traj_iniziale[1], traj_iniziale[2]), 2.) );
			sum_A += A;
	
			double t = rnd.Rannyu();
			if(t <= A){
				for(int n = 0; n < 3; ++n){
					traj_iniziale[n] = traj[n];
				}
			}

			sum_r += r_sferiche(traj_iniziale[0], traj_iniziale[1], traj_iniziale[2])/a_bohr;
		}
		r_medio[i] = sum_r / double(L);
		cout << "accettanza: " << sum_A/(double(L)) << endl;
	}

	fout2.close();
	MediaBlocchi(N, r_medio, "unif_r210", "unif_r210_err");

	return 0;
}
