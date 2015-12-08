#include "specie.h"

Specie::Specie()
{
	spname = "test";
}



void Specie::InitialData(int id, string sp_name, double* composition)
{
	ID = id;
	spname = sp_name;
	totalN = composition;

	//if(type == "atom")
	//{
	InitCon();
	InitLine();
	//}

}
double Specie::GetN(int T)
{
	if (T < compositionts)
		cout << "temperature must be greater than "<<ts<<" K" << endl<<"reference: "<<spname<<endl;
	 return totalN[(int)((T-compositionts)/compositiontstep)];
}

void Specie::AssignN(int T)
{
	N = GetN(T);
//	cout << N << endl;
}


void Specie::InitCon()
{
	ifstream datafile;  //file stream of data contains species
	string s;         //buffer string used to read string from txt.

	datafile.open(databasepath+"Continue\\"+spname+".txt",ios::in);
	if(datafile.fail())
		cout<<"open file failed"<<endl;

	datafile>>s;        //read in number

	if(isnum(s))
		n = atof(s.c_str());
	else
		cout<<"Wrong format in data file"<<endl;
	datafile>>s;       //read in net charge
	Za = atof(s.c_str());
	datafile>>s;      //read in ionization energy
	EI = atof(s.c_str());
	datafile>>s;
	Es = atof(s.c_str());
	J = new double[n];  //allocate J E g
	E = new double[n];
	g = new double[n];
	Neff = new double[n];
	int i = 0;
	nl = 0;
	while(!datafile.eof())
	{
		if(i==n)
			cout<<"out of range of n"<<endl;
		datafile>>s;
		if(isnum(s))
		{
			J[i] = atof(s.c_str());
			datafile>>s;
			E[i] = atof(s.c_str());  //read in data
			if(E[i]>=EI && nl==0)
				nl = i;
			g[i] = 2*J[i]+1;     //calculate while reading
			Neff[i] = sqrt(pow(Za,2)*EH/(EI-E[i])); 
			i++;
		}
		
	}
	if (nl == 0)
		nl = i;
	pf = 0;
	cout<<"Initial "<<spname<<" complete!"<<endl;
	datafile.close();
}

void Specie::InitLine()
{
	ifstream datafile;  //file stream of data contains species
	string s;         //buffer string used to read string from txt.
	
	datafile.open(databasepath + "Line\\" + spname + ".txt", ios::in);
	if(datafile.fail())
		cout<<"open file failed"<<endl;

	datafile>>s;        //read in number

	if(isnum(s))
		nline = atof(s.c_str());
	else
	{
		cout<<"Wrong format in data file"<<endl;
		return;
	}
	line = new class Line[nline];

	int i = 0;
	while(!datafile.eof())
	{
		if(i==nline)
			cout<<"out of range of n"<<endl;
		datafile>>s;
		if(isnum(s))
		{
			//line[i].seq = atof(s.c_str());
			//datafile>>s;
			line[i].fq = atof(s.c_str()); //line center frequency
			datafile>>s;
			line[i].wn = atof(s.c_str()); //wavenumber, energy of line
			datafile>>s;
			line[i].f = exp(atof(s.c_str()));  //read in data
			datafile>>s;
			line[i].le = atof(s.c_str());
			datafile>>s;
			line[i].lj = atof(s.c_str());
			datafile>>s;
			line[i].ue = atof(s.c_str());
			datafile>>s;
			line[i].uj = atof(s.c_str());
			line[i].wl = hkcm*c/line[i].wn;
			if (line[i].f > 0.1)
			i++;
		}
		
	}
	if (line[i].f <= 0.1)
		i = i - 1;
	nline = i;
	datafile.close();
}

void Specie::cal_nor(double vs,double vf, double vp, int linenum)
{
	
	
	int i = linenum;
		line[i].normalise = 0;
		for (double v = vs; v < vf; v += vp)
		{
			line[i].normalise += 1 / (PI*line[i].stark*(pow((v - line[i].fq) / line[i].stark, 2) + 1))*vp;
		}
		
	
	
	//line[i].normalise *= 2;
}
 
void Specie::cal_voigt_nor(double vs, double vf, double vp, int linenum)
{
	double sk = line[linenum].stark;
	line[linenum].voigtnormalise = 0;
	double fq = line[linenum].fq;
	
	for (double v = vs; v < vf; v += vp)
	{
		double X1 = (v - fq);
		G = exp(-pow(X1, 2) / (2 * pow(sigma, 2))) / (sigma*pow(2 * PI, 0.5));
		
		double L = 1 / (PI*sk*(pow(X1/sk, 2) + 1));
		line[linenum].voigtnormalise += (eta*L + (1 - eta)*G)*vp;
	}
		
}

void Specie::printallstark()
{
	ofstream startout(path + spname+"start.txt");
	for (int l = 0; l<nline; l++)
		startout << l << "  " << line[l].stark <<"   "<< line[l].normalise<<"  "<<line[l].fq<<endl;
	startout.close();
}

void Specie::cal_stark(int T,double Ne)
{
	double funcF = 16*pow((PI/3),2/3)*c*Ry*pow(Br,3)*Ne*pow(hk*c*Ry/(BolJ*T),0.5);
 	for(int l=0;l<nline;l++)
	{
		double s_up = 0;
		double s_low = 0;
		int i = 0;
		while (find_line(line[l].ue, line[l].uj,'u',i) != -1)
	{
		i = find_line(line[l].ue, line[l].uj, 'u', i);
		s_up += line[i].linestrength()*line[i].gse(T)/2*line[i].uj;
		//cout<<"UE"<<i<<"  "<<sp[n].line[i].ue<<"   "<<sp[n].line[i].uj<<endl;
		i++;
	}
	i = 0;
	while (find_line(line[l].le, line[l].lj, 'l', i) != -1)
	{
		i = find_line(line[l].le, line[l].lj, 'l', i);
		s_low += line[i].linestrength()*line[i].gse(T)/2*line[i].lj;
		//cout<<"LE"<<i<<"  "<<sp[n].line[i].le<<"   "<<sp[n].line[i].lj<<endl;
		i++;
	}
	line[l].stark = funcF*(s_up+s_low);
	if (line[l].stark == 0)
		cout << "test" << endl;
	/*if (line[l].stark > 1e12)
		cout << spname << "  " << l << " " << line[l].stark << endl;*/
	}
}

int Specie::find_line(double levelenergy, double levelj, char c, int startnumber)
{
	
	while(startnumber<nline)
	{
		
		for(int i=startnumber;i<nline;i++)
		{
			if (c == 'u' && line[i].ue == levelenergy &&line[i].uj == levelj)
			return i;
			if (c == 'l' && line[i].le == levelenergy &&line[i].lj == levelj)
			return i;
		}
		return -1;
	}
	return -1;

}

bool Specie::isnum(string s)  //pick only number
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



void Specie::cal_ln(int T)
{
	if(pf==0)
		calPF(T);
	double weight;
	for(int i=0;i<nline;i++)
	 {
		 weight = 2*line[i].lj+1;
		 line[i].ln = weight*N*exp(-(line[i].le/(Bolm*T)))/pf;
		// cout << line[i].ln << endl;
	}

}

void Specie::show()
{
	int i=0;
	cout<<"   J   "<<"  Energy  "<<endl;
	for(i;i<n;i++)
		cout<<J[i]<<"    "<<E[i]<<"   "<<g[i]<<endl;
	cout<<"n = "<<n<<endl<<"i = "<<i<<endl;
}

double Specie::calPF(int T) //calculate total pf
{
	pf = 0;
	for(int i = 0;i<n;i++)
		pf+=g[i]*exp(-E[i]/(Bolm*T));
	return pf;
}