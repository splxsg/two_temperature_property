#ifndef _SPECIE_H_
#define _SPECIE_H_

#include "parameter.h"
#include "Line.h"

class Specie : public Line
{
public:
	Specie();
	void InitialData(int id, string sp_name, double *composition);
	int ID;
	void InitialData();
	void ReadinComposition(double* comp);
	void show(); //for test display data
	string type;
	string spname;

	///continue part
	void InitCon();
	double Za;
	int n;
	int nl; //number less than limit;
	double pf; //total parittion function;
    double calPF(int T); //total partition function 
	double calquantumnumber;

	double GetN(int T);
	void AssignN(int T);
	double* totalN;
	double* J;
	double* E;
	double* g; //statistical weight
	double* Neff; //Effective principal quantum number
	double N; //number density
	
	
	double EI;
	double Es; //lowest smeared energy 
	bool isnum(string s); 
	

	//line part
	void InitLine();
	int nline;
	int find_line(double levelenergy, double levelj, char c, int startnumber);
	class Line* line;

	void cal_stark(int T,double Ne);
	void cal_ln(int T);
	//class Line line;
	void cal_nor(double vs, double vf,double vp, int linenum); //start, finish, step
	void printallstark();
	void cal_voigt_nor(double vs, double vf, double vp, int linenum);


};



#endif