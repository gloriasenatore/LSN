#ifndef __Funzioni__
#define __Funzioni__

void ToFile(const char* Filename, double * data, int size);
void ToFileApp(const char* Filename, double * data, int size);
double uncertainty(double ave2, double ave, int n);
void MediaBlocchi(int N, double * valori, double * valori2, const char * filename_valori, const char * filename_err);
void MediaBlocchi(int N, double * valori, const char * filename_valori, const char * filename_err);
double r_sferiche(double x, double y, double z);
double theta_sferiche(double x, double y, double z);
double phi_sferiche(double x, double y);
double psi_100(double x, double y, double z);
double psi_210(double x, double y, double z);

#endif
