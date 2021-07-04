#ifndef __Funzioni__
#define __Funzioni__

void ToFile(const char* Filename, double * data, int size);
void ToFileApp(const char* Filename, double * data, int size);
double uncertainty(double ave2, double ave, int n);
void MediaBlocchi(int N, double * valori, double * valori2, const char * filename_valori, const char * filename_err);
void MediaBlocchi(int N, double * valori, const char * filename_valori, const char * filename_err);
double psi_T(double sigma, double mu, double x);
double d2_psi_T(double sigma, double mu, double x);
double V(double x);

#endif
