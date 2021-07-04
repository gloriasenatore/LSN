/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  bool restart = false;     // set true per ripartire da una configurazione precedente
  bool equilibrate = false;  // set true per valutare l'equilibrazione del sistema
  Input(restart);           //Inizialization
  bool values = false;       // set true per stampare solo i risultati con medie progressive e incertezze finali, serve per plottare le osservabili in funzione della temperatura

  if(equilibrate == true){
     for(int istep=1; istep <= 1000; ++istep)
     {
       Move(metro);
       Measure();
       cout << istep << "\t" << "U = " << walker[iu] << "\t" << "M = " << walker[im] << endl;  
     }
  }

  for(int istep=1; istep <= 500; ++istep) Move(metro); // lascio equilibrare il sistema per 500 passi

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages

    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk, values);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(bool restart)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program performs Metropolis moves" << endl;
  else cout << "The program performs Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  iu2 = 4; //Energia interna al quadrato
  im2 = 5; //Magnetizzazione al quadrato
 
  n_props = 6; //Number of observables --> ho aggiunto due osservabili: energia interna al quadrato e magentizzazione al quadrato

//initial configuration
  if(restart == false){
     for (int i=0; i<nspin; ++i)
     {
       if(rnd.Rannyu() >= 0.5) s[i] = 1;
       else s[i] = -1;
     }
  }
//Start from a previous configuration
  if(restart == true){ cout << "restart" << endl;
       fstream Conf;
       Conf.open("config.final");

       for (int i=0; i<nspin; ++i){
          Conf >> s[i];
       }
       Conf.close();
  } 

  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm, acc;
  double p_up, p_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin); //Non prende tutti e 50 gli spin diversi

    if(metro==1) //Metropolis
    {
	energy_old = Boltzmann(s[o], o);
	p = s[o]*(-1);
	energy_new = Boltzmann(p, o);
	acc = min(1., exp(-beta*(energy_new - energy_old)));

	double r = rnd.Rannyu();
	if(r <= acc){
		s[o] = p;
		accepted++;
	}
	attempted++;
    }
    else //Gibbs sampling
    {
	p_up = exp(-beta*Boltzmann(1, o));
	p_down = exp(-beta*Boltzmann(-1, o));
	p = p_up / (p_up + p_down);

	double r = rnd.Rannyu();
	if(r <= p){
		s[o] = 1;
	} else s[o] = -1;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double ene = 0.0, u = 0.0, u2 = 0.0, m = 0.0, m2 = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i) //sommo su tutti gli spin
  {
     ene = -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]); //Energia interna
     u += ene;
     u2 += pow(ene, 2.); //Energia interna al quadrato
     m += s[i]; //Magnetizzazione
     m2 += pow(s[i], 2.); //Magnetizzazione al quadrato
  }
  walker[iu] = u;
  walker[iu2] = u2;
  walker[im] = m;
  walker[im2] = m2;
}


void Reset(int iblk) //Reset block averages
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
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i) //blk_av è lungo m_props = 1000 perché possiamo misurare fino a 1000 osservabili
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk, bool values) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi, Values;
   const int wd=15;
    
    cout << "Block number " << iblk << endl;
    if(metro == 1) cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    if(values == false) Ene.open("output.ene."+to_string(temp),ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u; //Media progressiva
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    if(values == false) Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    if(values == false) Ene.close();

    if(values == false) Heat.open("output.heat."+to_string(temp),ios::app);
    stima_c = pow(beta, 2.) * (blk_av[iu2]/blk_norm/(double)nspin - pow(blk_av[iu]/blk_norm/(double)nspin, 2.)); //Capacità termica
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    if(values == false) Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    if(values == false) Heat.close();

    if(values == false) Mag.open("output.mag."+to_string(temp),ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetizzazione
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    if(values == false) Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    if(values == false) Mag.close();  

    if(values == false) Chi.open("output.chi."+to_string(temp),ios::app);
    stima_x = beta * blk_av[im2]/blk_norm/(double)nspin;; //Suscettività
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    if(values == false) Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    if(values == false) Chi.close();  

    if(iblk == nblk){ // Stampa i risultati considerando tutti i blocchi
    	Values.open("values.dat",ios::app);
    	Values << setw(wd) << temp << setw(wd) << "ene" << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << setw(wd) << "heat" << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << setw(wd)
	<< "mag" << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << setw(wd) << "chi" << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    	Values.close();
    }       

    cout << "----------------------------" << endl << endl;
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
