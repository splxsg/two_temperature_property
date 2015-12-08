#include "speciem.h"

Speciem::Speciem()
{

}

void Speciem::InitialData(string sp_name)
{
	InitCS(sp_name);

}

bool Speciem::isnum(string s)  //pick only number
{
	stringstream sin(s);
	double t;
	char p;
	if(!(sin >> t))
		return false;
	if(sin >> p)
		return false;
	else
		return true;
}


void Speciem::InitCS(string sp_name)
{
	ifstream datafile;  //file stream of data contains species
	string s;         //buffer string used to read string from txt.

	datafile.open(path+"continue\\M\\"+sp_name+".txt",ios::in);
	if(datafile.fail())
		cout<<"open file failed"<<endl;

	datafile>>s;        //read in number

	if(isnum(s))
		v12 = atof(s.c_str());
	else
		cout<<"Wrong format in data file"<<endl;
	datafile>>s;       //23 frequency
	v23 = atof(s.c_str());
}

double Speciem::GetN(int T)
{
	return totalN[(int)((T-2000)/10)];
}

void Speciem::AssignN(int T)
{
	N = GetN(T);
}
