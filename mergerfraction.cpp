#include "halomass.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>
#include <time.h>
#include <iomanip>
//#include <TF1.h>
//#include "brent.hpp"
//using namespace brent;

#include <boost/multi_array.hpp>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/normal_distribution.hpp>
//#include <boost/random/variate_generator.hpp>

using namespace std;
typedef vector<double>::iterator iter;

struct LambdaCDM
{
	//Parameters for a Lambda-CDM cosmology with zero curvature from Planck 2013 results (http://arxiv.org/pdf/1303.5076v1.pdf)
	//Changed to WMAP 7 cosmology instead to be in accordance with MPA-JHU data
	LambdaCDM(): H_0(70.0), omega_M(0.3), omega_l(0.7) {}
	double H_0;
	double omega_M;
	double omega_l;
} LCDM;

bool GetMergerCoefficients(double (&x)[5], string);
void RunModel(double, double x[]);


/*double RedshiftToTime(double redshift)
{
	//Solve t=int(dz)_0^z by simple numerical integration (time in Gyr), z(0) defined to be approx. 13.811 Gyr
	double delta_z = 0.000001;
	double t = 0;
	for (double z=0;z <= redshift;z += delta_z)
	{
		t += delta_z/((1+z)*sqrt(LCDM.omega_M*pow(1+z, 3)+LCDM.omega_l));
	}
	t *= (1/((LCDM.H_0 * 1000/3.08567758e22))/3.15576e16);
	return t;
}

double RedshiftToTimeZero(double redshift)
{
	//Zero function needed for root finding algorithm
	return RedshiftToTime(redshift) - curtime;
}

double TimeToRedshift2(double time)
{

	//Use Brent root finding algorithm to invert t=int(dz)_0^z (time in Gyr), t(13.811) defined to be 0
	double bracketl = 0;
	double bracketh = 10;
	double tolerance = r8_epsilon();
	curtime = time;
	return zero(bracketl, bracketh, tolerance, RedshiftToTimeZero);
}*/


class Mergers
{
	double x_min;
	double t_0;
	double TimeToRedshift(double);
	double dMergers(double, double, double);

public:
	Mergers(double x, double t): x_min(x), t_0(t) {}

	double NMergers(double, double, double, double);
	double dMergersTimeIntegrated(double, double,  double);
};

double Mergers::TimeToRedshift(double time)
{

	//Use Brent root finding algorithm to invert t=int(dz)_0^z (time in Gyr), t(13.811) defined to be 0
	/*double bracketl = 0;
	double bracketh = 10;
	double tolerance = r8_epsilon();
	curtime = time;
	return zero(bracketl, bracketh, tolerance, RedshiftToTimeZero);*/

	//Use low-redshift linear approximation for speed
	return 0.071*time+0.0045*time*time;//for Planck cosmology: 0.069*time+0.0006*time*time;
}

//For halo merger rate using halo mass ratio
double Mergers::dMergers(double massratio, double HaloMassTerm, double time)
{
	//Calculate galaxy merger rate according to Fakhouri et al. 2010
	static double InverseHubbleTimeGyr = (LCDM.H_0 * 1000/3.08567758e22)*3.15576e16; //Hubble parameter in units of Gyr^-1
	double z = TimeToRedshift(time);
	//double dNdxdz = 0.0104 * pow(HaloMass / 1e12, 0.133) * pow(massratio, -1.995)*exp(pow(massratio/9.72e-3,0.263)) * pow(1+z, 0.0993);
	double dNdxdz = 0.0104 * HaloMassTerm * pow(massratio, -1.995)*exp(pow(massratio/9.72e-3,0.263));
	//Convert dN/dxdz to dN/dxdt by dN/dxdt=dN/dxdz * dz/dt and dz/dt = H(z)*(1+z) (frequency in Gyr^-1)
	//(1+0.51*z) is an approximation of sqrt(LCDM.omega_M*pow(1+z, 3)+LCDM.omega_l) for omega_M=0.315 and omega_l=0.685
	//(1+0.48*z) is an approximation of sqrt(LCDM.omega_M*pow(1+z, 3)+LCDM.omega_l) for omega_M=0.3 and omega_l=0.7
	//These are valid for z < 0.2
	//If generality is prefered to speed of calculation replace this by sqrt(LCDM.omega_M*pow(1+z, 3)+LCDM.omega_l)
	double dNdxdt = dNdxdz * (1+z) * InverseHubbleTimeGyr * (1+0.48*z);//(1+0.51*z);
	return dNdxdt;
}

double Mergers::dMergersTimeIntegrated(double Delta_t, double HaloMassTerm,  double x)
{
	double t_span = 3.42;		//Span of time in Gyr, approximate cosmic time at z=0.3
	double dNmergerdx = 0;

	for (double t = 0; t <= t_span; t += Delta_t) //loop over "time before present" bins (time in Gyr)
	{
		//double M_hs = M_hb * x;
		double dNmergerdxdt = dMergers(x, HaloMassTerm, t_0 - (t - t_span));
		dNmergerdx += Delta_t * dNmergerdxdt;
	}
	
	return dNmergerdx;
}

double Mergers::NMergers(double Delta_t, double x_min, double Delta_x, double HaloMassTerm)
{
	double R = 0;
	for(double x = x_min; x + Delta_x <= 1; x += Delta_x)
	{
		double R1 = dMergersTimeIntegrated(Delta_t, HaloMassTerm, x);
		double R2 = dMergersTimeIntegrated(Delta_t, HaloMassTerm, x + Delta_x);
		R += 0.5 * (R1 + R2) * Delta_x;
	}
	return R;
}

int main(int argc, char* argv[])
{
	string FMRtype = "O3N2MarinoFMR-Brinchmann";
	double MergerWeights[5];
	double x_min = 0.205;
	if (!GetMergerCoefficients(MergerWeights, FMRtype)) {cout << "Failed to load coeffiecients for Gaussian - Quitting." << endl;return -1;}

	RunModel(x_min, MergerWeights);

	return 0;
}


//Read merger mass weights from file
bool GetMergerCoefficients(double (&Weights)[5], string Type)
{
	ifstream datain("FMRCoeffs.txt");
	if (!datain.is_open()) {cout << "Unable to open file FMRCoeffs.txt"<< endl;return false;}
	string line;
	int i = 0;
	bool typefound = false;
	bool gaussfound = false;
	while (getline(datain, line))
	{
		if (line == "Type: " + Type) {typefound=true;continue;}
		if (!typefound) {continue;}
		if ((typefound) && (line.substr(0,5) == "Gauss")) {gaussfound = true;continue;}
		if (gaussfound)
		{
			if (i >= 3) {Weights[i-3] = atof(line.c_str());}
			if (i == 2) {getline(datain, line);i++;continue;}
			if (i == 7) {break;}
			i++;
		}
	}
	datain.close();
	if (i != 7) {return false;} else {return true;}
}

void RunModel(double x_min, double MergerWeights[])
{
	//const double t_now = 13.811;
	double Delta_t = 0.001;
	//double Delta_x = 0.01;
	const double Delta_x_fine = 1e-4;
	//const double NMergerTolerance = 1e-5;
	//double t_0 = 0;		//Starting time in Gyr
	//double M_sb = 1e12;	//Stellar mass of big galaxy
	Mergers Merger(x_min, 0.93);
	double MergerFraction = 0;

	int i = 0;
	for (double M_sb = 9.125;M_sb < 11.5; M_sb += 0.25)
	{
		cout << M_sb << endl;
		double M_hb = HaloMass(pow(10, M_sb));	//Halo mass of big galaxy
		//Delta_x = 0.000390625;
		double HaloMassTerm = pow(M_hb / 1e12, 0.133);
		MergerFraction += Merger.NMergers(Delta_t, x_min, Delta_x_fine, HaloMassTerm) * MergerWeights[i];
		i++;
		//cout << "thread " << omp_get_thread_num() << " T" << i+1 << ": " << omp_get_wtime()-t_start << endl;
	}

	cout << "Merger fraction: " << MergerFraction << endl;
}
