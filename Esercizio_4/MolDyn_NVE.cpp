/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  bool restart = true;         // set restart false se è la prima partenza, altrimenti true se si sta ripartendo
  bool equilibrate = false;     // set equilibrate true se si vuole equilibrare il sistema
  bool blocking = true;        // set blocking true se si sta operando il blocking average, altrimenti false. L'equilibrazione avviene solo se blocking == false e equilibrate == true.
  Input(restart);              //Inizialization
  int nconf = 1;

  if(blocking == false){

   for(int istep=1; istep <= nstep; ++istep){

    if(equilibrate == true){
      if (istep == 1 || istep == 1000 || istep == 3000 || istep == 5000 || istep == 7000) Equilibration(); //scelgo di riscalare le velocità per cinque step separati da passi dove il sistema si muove normalmente con Verlet
      else Move();           //Move particles with Verlet algorithm
   }

    else Move();            //Move particles with Verlet algorithm

     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure(blocking);     //Properties measurement
//      ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
   }
  }

  else{

      for(int iblk=1; iblk <= 100; ++iblk) //scelgo un numero di blocchi pari a 100
      {
         Reset(iblk);

         for(int istep=1; istep <= 1000; ++istep) //ci sono 1000 in ogni blocco
         {
         	Move();		
		if(istep%10 == 0) Measure(blocking);		
      		Accumulate(); 
         }
         Averages(iblk);
       }
  }

  ConfFinal();         //Write final configuration to restart
  ConfFinal_1(); //Write configuration t-delta t to restart

  return 0;
}


void Input(bool restart){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  double pi = M_PI;

  //Tail corrections for potential energy and pressure
  //vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  //cout << "Tail correction for the potential energy = " << vtail << endl;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   if(restart == false){
   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   	double sumv[3] = {0.0, 0.0, 0.0};
   	for (int i=0; i<npart; ++i){
     		vx[i] = rand()/double(RAND_MAX) - 0.5;
     		vy[i] = rand()/double(RAND_MAX) - 0.5;
     		vz[i] = rand()/double(RAND_MAX) - 0.5;

     		sumv[0] += vx[i];	
     		sumv[1] += vy[i];
     		sumv[2] += vz[i];
   	}
   	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   	double sumv2 = 0.0, fs;
   	for (int i=0; i<npart; ++i){
     		vx[i] = vx[i] - sumv[0];
     		vy[i] = vy[i] - sumv[1];
     		vz[i] = vz[i] - sumv[2];

     		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   	}
   	sumv2 /= (double)npart;

   	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   	for (int i=0; i<npart; ++i){
     		vx[i] *= fs;
     		vy[i] *= fs;
     		vz[i] *= fs;

     		xold[i] = Pbc(x[i] - vx[i] * delta);
     		yold[i] = Pbc(y[i] - vy[i] * delta);
     		zold[i] = Pbc(z[i] - vz[i] * delta);
   	}
   }

// ESERCIZIO 4.1
   if(restart == true){
	fstream Conf;
  	Conf.open("config.final-1");

  	for (int i=0; i<npart; ++i){
		Conf >> xold[i] >> yold[i] >> zold[i];
		xold[i] = xold[i] * box;
    		yold[i] = yold[i] * box;
    		zold[i] = zold[i] * box;
  	}
  	Conf.close();
   } 

   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(bool blocking){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  if(blocking == false){ //se non si sta eseguendo la media a blocchi, scrive i valori di una sola simulazione
     Epot.open("output_epot.dat",ios::app);
     Ekin.open("output_ekin.dat",ios::app);
     Temp.open("output_temp.dat",ios::app);
     Etot.open("output_etot.dat",ios::app);
  }

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     
     //g(r)
     if(dr <= box/2.){	
	for(int n = 0; n < 100; ++n){
		if(n*bin_size < dr && dr <= (n+1)*bin_size){
			bin = n;
			walker[igofr+bin] = walker[igofr+bin] + 2.;
		}
	}
     }
     

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    walker[iv] = v; //aggiorno i walker per le calcolare i valori medi con il blocking average
    walker[ik] = t;
    walker[ie] = t+v;
    walker[it] = (2.0 / 3.0) * t/(double)npart;

    if(blocking == false){
       Epot << stima_pot  << endl;
       Ekin << stima_kin  << endl;
       Temp << stima_temp << endl;
       Etot << stima_etot << endl;

       Epot.close();
       Ekin.close();
       Temp.close();
       Etot.close();
    }

    return;
}

// 4.1 EQUILIBRATION
void Equilibration(void){ // Equilibrate with Verlet algorithm
  double xnew[m_part], ynew[m_part], znew[m_part], fx[m_part], fy[m_part], fz[m_part], vx2[m_part], vy2[m_part], vz2[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  double t = 0.;
  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew[i] = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew[i] = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew[i] = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew[i] - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew[i] - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew[i] - zold[i])/(2.0 * delta);

    vx2[i] = Pbc(xnew[i] - x[i])/delta; //Velocità al tempo t + dt/2
    vy2[i] = Pbc(ynew[i] - y[i])/delta;
    vz2[i] = Pbc(znew[i] - z[i])/delta;

    t += vx2[i]*vx2[i] + vy2[i]*vy2[i] + vz2[i]*vz2[i]; // Temperatura al tempo t + dt/2
    }

  t /= (double)npart;
  for(int i=0; i<npart; ++i){

    vx[i] = vx[i] * sqrt(3. * temp/t); // Riscalo le velocità al tempo t
    vy[i] = vy[i] * sqrt(3. * temp/t);
    vz[i] = vz[i] * sqrt(3. * temp/t);

    x[i] = Pbc(xnew[i] - delta * vx[i]); // Aggiorno le posizioni al tempo t
    y[i] = Pbc(ynew[i] - delta * vy[i]);
    z[i] = Pbc(znew[i] - delta * vz[i]);

    xold[i] = x[i]; // Uso le posizioni al tempo t aggiornate e quelle al tempo t + dt per ricominciare la simulazione
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew[i];
    y[i] = ynew[i];
    z[i] = znew[i];
  }
  return;
}

//4.2 BLOCKING
void Reset(int iblk) //Azzero i vettori contenenti le somme dei valori istantanei e progressivi delle osservabili termodinamiche
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
}

void Accumulate(void) //Aggiorno i valori delle osservabili termodinamiche
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Scrivo i valori progressivi di medie e incertezze delle osservabili termodinamiche
{
   double r; 
   ofstream epot, ekin, etot, temp, gave, pres;
   const int wd=15;
    
    cout << "Block number " << iblk << endl;
    
    epot.open("ave_epot.out",ios::app);
    stima_ep = blk_av[iv]/blk_norm/(double)npart /*+ vtail*/; //energia potenziale per particella
    glob_av[iv]  += stima_ep; //media progressiva
    glob_av2[iv] += stima_ep*stima_ep;
    err_ep = uncertainty(glob_av2[iv]/(double)iblk, glob_av[iv]/(double)iblk, iblk-1);
    epot << setw(wd) << iblk <<  setw(wd) << stima_ep << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_ep << endl;
    epot.close();

    ekin.open("ave_ekin.out",ios::app);
    stima_ek = blk_av[ik]/blk_norm/(double)npart; //energia cinetica per particella
    glob_av[ik]  += stima_ek; //media progressiva
    glob_av2[ik] += stima_ek*stima_ek;
    err_ek = uncertainty(glob_av2[ik]/(double)iblk, glob_av[ik]/(double)iblk, iblk-1);
    ekin << setw(wd) << iblk <<  setw(wd) << stima_ek << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_ek << endl;
    ekin.close();

    etot.open("ave_etot.out",ios::app);
    stima_et = blk_av[ie]/blk_norm/(double)npart; //energia totale per particella
    glob_av[ie]  += stima_et; //media progressiva
    glob_av2[ie] += stima_et*stima_et;
    err_et = uncertainty(glob_av2[ie]/(double)iblk, glob_av[ie]/(double)iblk, iblk-1);
    etot << setw(wd) << iblk <<  setw(wd) << stima_et << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_et << endl;
    etot.close();

    temp.open("ave_temp.out",ios::app);
    stima_t = blk_av[it]/blk_norm; //temperatura
    glob_av[it]  += stima_t; //media progressiva
    glob_av2[it] += stima_t*stima_t;
    err_t = uncertainty(glob_av2[it]/(double)iblk, glob_av[it]/(double)iblk, iblk-1);
    temp << setw(wd) << iblk <<  setw(wd) << stima_t << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_t << endl;
    temp.close();

    //confronto con esercitazione 7
    gave.open("ave_g.out",ios::app);
    for (int k=igofr; k<igofr+nbins; ++k){
	r = (double)(k-4) * bin_size;
	stima_g = blk_av[k]/blk_norm/npart/rho/(4.*M_PI/3.)/(pow(r + bin_size, 3.) - pow(r,3.));
	glob_av[k] += stima_g;
	glob_av2[k] += stima_g*stima_g;
	err_gdir = uncertainty(glob_av2[k]/(double)iblk, glob_av[k]/(double)iblk, iblk-1);

        if(iblk == 100){
		gave << setw(wd) << r+bin_size/2. << setw(20) << glob_av[k]/(double)iblk << setw(20) << err_gdir << endl;
	}
    }
    gave.close();

    cout << "----------------------------" << endl << endl;
}

double uncertainty(double ave2, double ave, int n){

	if(n == 0){
		return 0;
	}

	else{
		return sqrt((ave2 - pow(ave, 2)) / double(n));
	}
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfFinal_1(void){ //Write configuration at time t - delta t
  ofstream WriteConf;

  cout << "Print configuration at time t - delta t to file config.final-1 " << endl << endl;
  WriteConf.open("config.final-1");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
