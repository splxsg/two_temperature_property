#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#include "iostream"
#include "string"
#include "fstream"
#include <math.h>
#include <sstream>


using namespace std;

extern string path; //path of database
extern string pressure;
extern int p1; //first specie %
extern int p2; //second specie %
extern int p3; //third specie %
extern double BolJ; // Boltzmann constant J/K
extern double Bolm; //  Boltzmann cm-1/K
extern double Br;  //  Bohr radius m
extern double al; // Fine-structure constant
extern double hk; //Planck constant  m2kg/s
extern double PI; // PI
extern double c; // light speed
extern double EH; //Ionization energy of hydrogen  cm-1
extern double hkcm; //Planck constant cm-1.s
extern double Ry; //Rydberg constant m-1
extern double Ro; //classical electron radius
extern double C1[5]; //polyfit coefficient for gse;
extern double C2[5]; //polyfit coefficient for gsh;
extern string databasepath;
extern string path;
extern double ts, tf, tstep,compositionts,compositiontf,compositiontstep;
extern double sigma, gamma, FL, FG, f1, eta, G;
#endif