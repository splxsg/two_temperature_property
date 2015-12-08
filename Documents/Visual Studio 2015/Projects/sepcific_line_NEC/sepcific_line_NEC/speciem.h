#ifndef _SPECIEM_H_
#define _SPECIEM_H_

#include "parameter.h"

class Speciem
{
public:
	Speciem();
	void InitialData(string sp_name);
	void InitCS(string sp_name);


	//cross section
	double v12; //threshold frequency of region 1-2 
	double v23; //threshold frequency of region 2-3

	bool isnum(string s);
	double GetN(int T);    //get number density of specific temperature
	void AssignN(int T);     //assign number density to N, 
	double* totalN;             //array with all number density
	double N; //number density, must assign new number density when T changed
	
};

#endif