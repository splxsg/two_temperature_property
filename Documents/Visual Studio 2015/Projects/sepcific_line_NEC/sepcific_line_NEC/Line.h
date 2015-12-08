#ifndef _LINE_H_
#define _LINE_H_

#include "parameter.h"

class Line
{
public:
		Line();
		int seq; //NO. of line
	double wl; //wavelength
	double fq; //line center frequency
	double wn; //wavenumber
	double f; //oscillator strength
	double le; //lower level energy
	double lj; //lower level J
	double ue; //upper level energy
	double uj; //upper level J
	double ln; //number density of lower band
	double linestrength();
	double gse(int T);
	double gsh(int T);
	double normalise;
	double stark;
	double voigtnormalise;

};
#endif