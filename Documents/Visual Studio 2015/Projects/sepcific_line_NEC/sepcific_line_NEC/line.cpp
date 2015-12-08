#include "line.h"

Line::Line()
{
	normalise = 0;
	stark = 0;
	voigtnormalise = 0;
}

double Line::linestrength()
{
	return 3*(2*lj+1)*wl*Ry*f;
}

double Line::gse(int T)
{
	double x = Bolm*T*3/(2*(wn));
	return pow(x,4)*C1[0]+pow(x,3)*C1[1]+pow(x,2)*C1[2]+pow(x,1)*C1[3]+C1[4];
}

double Line::gsh(int T)
{
	double x = Bolm*T*3/(2*(wn));
	if(x>=20)
		return 1e-3*(x-20)+0.81;
	else
		return pow(x,4)*C2[0]+pow(x,3)*C2[1]+pow(x,2)*C2[2]+pow(x,1)*C2[3]+C2[4];
}
