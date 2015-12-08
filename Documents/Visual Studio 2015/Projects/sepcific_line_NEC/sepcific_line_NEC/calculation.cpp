#pragma once
#include "calculation.h"





int main()
{
	double ni, cs, Kc, vp, Ksm, C1, C2, C3, C4;
	
	//ask for profile number
	string pathst;
	cin >> pathst; //if multiple tasks required, uncomment this and comment next one
	//pathst = "11";
	path += "case" + pathst +"\\";
	//ask for profile number

	//assign output file
	ofstream fsc(path + "abs.txt");
	
	ofstream fout(path + "NEC.txt");
	//assign output file

	SystemInit();
	/*for (int s = 0; s < N_atom; s++)
		cout << sp[s].spname<<"  "<<sp[s].nl << endl;*/
	
	


		Kci = new double[N_sp];
		Ksm = 0;
		double Ne;
		double Kl = 0;
		double templine = 0;
		double voigtline = 0;
		double vs = 2e12;
		double vf = 1e16;
		double step = 1e12;
		int linenumber = 21;
		int speicenumber = 10;
	
		
		
			
		cout << "Start cal" << endl;
		
				
			for (T = ts; T < tf; T += tstep){
				Ne = tNe[(int)((T - compositionts) / compositiontstep)];
				double *absK = new double[int((vf - vs) / step)];

				
					sp[speicenumber].AssignN(T);
					sp[speicenumber].cal_stark(T, Ne);
					sp[speicenumber].cal_ln(T);
					sp[speicenumber].cal_nor(vs, vf, step, linenumber);
					
					//sp[i].printallstark();

				
			//	cout << sp[0].line[0].fq << endl;
				cout << "============================================================================================" << endl;
				cout << "Temperature is " << T <<endl;
				cout << "Progress: ";
				string absfile = path + to_string(T) + ".txt";
				string absconfile = path + to_string(T) + "_voigt.txt";
				ofstream abskfout(absfile);
				ofstream abskconfout(absconfile);

///test
				int s = speicenumber;
				int i = linenumber;
				double templinesum = 0;
				double voigtsum = 0;
				sigma = sp[s].line[i].stark / 1.1;
				gamma = sp[s].line[i].stark;
				FL = 2.0 * gamma;
				FG = 2.0 * sigma*sqrt(2 * log(2));
				f1 = pow(FG, 5) + 2.69269*pow(FG, 4)*FL + 2.42843*pow(FG, 3)*pow(FL, 2) + 4.47163*pow(FG, 2)*pow(FL, 3) + 0.07842*FG*pow(FL, 4) + pow(FL, 5);
				f1 = pow(f1, 0.2);
				eta = 1.36603*(FL / f1) - 0.47719*pow(FL / f1, 2) + 0.11116*pow(FL / f1, 3);
				G = 0;
				sp[speicenumber].cal_voigt_nor(vs, vf, step, linenumber);
				for (double v = vs; v < vf; v += step)
				{
					

					
					Kl = 0;
					Kc = 0;
					C1 = (BolJ*T*pow(al, 3)*pow(c, 2)) / (4 * hk*pow(v, 3));
					C4 = 0;
					templine = 0;
	
					

					
					//for (int s = 1; s < 2; s++)
					
						
						Kci[s] = 0;
						double vn = (sp[s].EI - sp[s].Es) / hkcm;
						C2 = sp[s].N / sp[s].pf*pow(sp[s].Za, 2);
						//continue part
						bool flag = true;
						C3 = 0;
						
						
							
							if ( v >= (sp[s].EI - sp[s].E[i]) / hkcm && flag)
							{
								ni = NI(sp[s].N, sp[s].g[i], sp[s].E[i], T, sp[s].pf);
								//cout << sp[s].Neff[i] << endl;
								cs = CS(sp[s].Neff[i], (sp[s].EI - sp[s].E[i]) / hkcm, sp[s].Za, v);
								Kci[s] += ni*cs;
								//Kci[s] += cs;
								flag = false;
							}

							double sk = sp[s].line[i].stark;
								
							templine = LorentzP(v, T, s, i, sk);
				
							double X = (v - sp[s].line[i].fq);
							G = exp(-pow(X, 2) / (2 * pow(sigma, 2))) / (sigma*pow(2 * PI, 0.5));
							voigtline = voigtP(templine, G, eta)/sp[s].line[i].voigtnormalise;
							templine /= sp[s].line[i].normalise;
							templinesum += templine;
							voigtsum += voigtline;
							

									
								if (Kl < templine)
							Kl = templine;
							if (double::IsNaN(Kl))
								cout << "test" << endl;
							
						
					

					if (fmod(v, 1e15) == 0.0)
						cout << v << "    ";
			
					
					
					abskfout << v << "  " <<Kl << endl;
					abskconfout << v << "  " << voigtline<< endl;
				}

				abskfout.close();
				abskconfout.close();

				
				
				
				fout << T << "   ";
				
				
				cout << endl;
				cout << "templinesum = " << templinesum << endl;
				cout << "normalise = " << sp[s].line[i].normalise << endl;
				cout << "viogtsum = " << voigtsum << endl;
				cout << "normalise = " << sp[s].line[i].voigtnormalise << endl;
				cout << "=============================================================================================" << endl << endl;
				fout << endl;
			}
	
			fout.close();
				

	system("pause");
}
double f_nec(double nec, double T)
{
	double coef[8] = { 4.58990614012818e-27, -5.83593481533401e-22, 3.14468158905234e-17, -9.29858864822877e-13, 1.62762973088971e-08, -0.000168493932149844, 0.954790318421663, -2285.02391071045 };
	nec *= 0.4;
	//cout << ra << "-" << nec << " ";
	if (T >= 12000 && T <= 24000)
		nec -= 1e11*(pow(T, 7)*coef[0] + pow(T, 6)*coef[1] + pow(T, 5)*coef[2] + pow(T, 4)*coef[3] + pow(T, 3)*coef[4] + pow(T, 2)*coef[5] + pow(T, 1)*coef[6] + pow(T, 0)*coef[7]);
	if (T >= 23000)
		nec *= exp(-(T - 23000.0) / 9100.0);
	return nec;
}
void split(string *dest, string str, int n, char *delims)
{
	char *result = NULL;
	char *cstr = new char[str.length() + 1];
	int i = 1;
	strcpy(cstr, str.c_str());
	result = strtok(cstr, delims);
	dest[0] = result;
	while (i < n && result != NULL)
	{
		result = strtok(NULL, delims);
		dest[i++] = result;
	}
}



void loadparameter()
{
	string spgroup;
	config opc(path + "config.ini");
	string* refstr;
	N_sp = atoi(opc.getValue("spn").c_str());
	ts = atoi(opc.getValue("nects").c_str());
	tf = atoi(opc.getValue("nectf").c_str());
	tstep = atoi(opc.getValue("nectstep").c_str());
	compositionts = atoi(opc.getValue("compositionts").c_str());
	compositiontf = atoi(opc.getValue("compositiontf").c_str());
	compositiontstep = atoi(opc.getValue("compositiontstep").c_str());


	N_T = (int)((compositiontf - compositionts) / (compositiontstep)) + 1;
	totalspe = new string[N_sp];
	spetype = new string[N_sp];
	//cout <<" test" <<opc.getValue("type") << endl;
	spcomposition = new double*[N_T];
	for (int i = 0; i < N_sp; i++)
		spcomposition[i] = new double[N_T];
	int temp;
	//ifstream composition_fin(path + "composition_formated.txt");
	ifstream composition_fin(path + "composition.txt");
	for (int i = 0; i < N_T; i++)
		for (int j = 0; j <= N_sp; j++)
			if (j == 0)
				composition_fin >> temp;
			else
			{
				composition_fin >> spcomposition[j - 1][i];
			//	cout << spcomposition[j - 1][i] <<"  ";
			}
	composition_fin.close();
	split(totalspe, opc.getValue("sp"), N_sp, ",[]");
	split(spetype, opc.getValue("type"), N_sp, ",[]");
	

	for (int i = 0; i < N_sp; i++)
		if (spetype[i] == "a")
			N_atom++;
	for (int i = 0; i < N_sp; i++)
		if (spetype[i] == "m")
			N_mole++;

	
	tNe = new double[N_T];
	
	

}

void SystemInit()
{
	loadparameter();
	sp = new class Specie[N_atom];
	spm = new class Speciem[N_mole];

	


	for (int i = 0, j = 0; i < N_sp; i++)
	{
		if (spetype[i] == "a")
			sp[j++].InitialData(i, totalspe[i], spcomposition[i]);
	}
	
	
	
	/*for (int i = 0, j = 0; i < N_sp; i++)
	{
		if (spetype[i] == "m")
			spm[i].InitialData(i, spe[j++], spcomposition[i + 1]);
	}*/
			
	
	
	for (int i = 0; i < N_sp; i++)
	{
		if (spetype[i] == "e")
			for (int j = 0; j < N_T; j++)
				tNe[j] = spcomposition[i][j];
	}
	cout << "e- Complete" << endl;

	cout << "All initial complete!" << endl;

}


double LorentzPC(double f,double n)
{
	//cout<<PI*Ro*c*f*n<<endl;
	return PI*Ro*c*f*n;
}
double voigtP(double L, double G, double eta)
{
	return eta*L + (1 - eta)*G;
}
double LorentzP(double v, double T, int n,int i,double sk)
{
	double fq = sp[n].line[i].fq;
	double X = (v-fq)/sk;
	//cout<<sp[n].line[i].fq<<endl;
	return 1/(PI*sk*(pow(X,2)+1));
}




//double Resonance()
double CS(double ni, double freqi, double Za, double freq)
{
	return al*pow(Br,2)*4*PI*PI*ni*ni/pow(Za,2)*pow(freqi/freq,3); //al*a^2*4*pi^2*Ni/Za^2*(Vi/V)^3
}

double NI(double N,double g,double E,double T, double pf)
{
	return g*N*exp(-E/(Bolm*T))/pf;
}


//void Initcomposition(double* tNe)
//{
//	int a = 3290;
//	stringstream ss1, ss2;
//	ss1 << p1;
//	string f1 = ss1.str();
//	ss2 << p2;
//	string f2 = ss2.str();
//	ifstream datafile;
//	datafile.open(path + "composition\\" + pressure + "\\" + f1 + "_" + f2 + "\\e-.txt", ios::in);
//	for (int i = 0; i < a; i++)
//		datafile >> tNe[i];
//}




