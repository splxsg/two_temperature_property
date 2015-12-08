#ifndef _CALCULATION_H_
#define _CALCULATION_H_

#include "specie.h"
#include "speciem.h"
#include "config.h"

class Specie* sp;
class Speciem* spm;
string databasepath = "D:\\Dropbox\\my document\\NEC2015\\spectra_database\\";
string path = "D:\\Dropbox\\my document\\NEC2015\\NEC\\N2PTFE\\";
string pressure = "1bar";
double compositionts;
double compositiontf;
double compositiontstep;
double CS(double N, double freqi, double Za, double freq);
double NI(double N,double g,double E,double T, double pf);
double Stark(double T,double Ne,int n, int l);
double LorentzP(double v, double T, int n,int i,double sk);
double LorentzPC(double f,double n);
double voigtP(double L, double G, double eta);
//void Initcomposition(double* tNe);
void loadparameter();
void split(string *dest, string str, int n, char *delims);
void SystemInit();

double f_nec(double nec, double T);
int N_sp, N_atom, N_mole, N_T;
double **spcomposition, *tNe;
string *totalspe, *spetype, *spe, *spem;
double *Kci;
double ts, tstep, tf;
double sigma, gamma, FL, FG, f1, eta, G;
int T;

#endif