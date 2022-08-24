//Calculate halo mass given stellar mass using equation 21 from BEHROOZI, CONROY & WECHSLER 2010 (BWC10)
//and stellar mass given halo mass using BEHROOZI, WECHSLER & CONROY 2013 (BCW 2013)

#ifndef HALOMASS_HPP_
#define HALOMASS_HPP_

#include <cmath>

double StellarMassFitFunc(double);

struct HaloParams
{
	//Default values with error estimates from BWC10 table 2 (at z=0), see that paper for interpretations of the parameters
	double beta;	//0.44 +0.04 -0.06
	double delta;	//0.57 +0.15 -0.06
	double gamm;	//1.56 +0.12 -0.38
	double M_1;		//12.35 +0.07 -0.16, equivalent to log10 of M_1 as it is given in eq 21 of BCW10 but equal to value given in table 2 for M_1 (confusing)
	double M_ss0;	//log10(M_ss0)=10.72 +0.22 -0.29, equivalent to log10 of M_s0 as it is given in eq 21 of BCW10 but equal to value given in table 2 for M_s0 (confusing)

	//Default values with error estimates from BCW13 section 5 (at z=0, a=1), see that paper for interpretations of the parameters
	//gamma, delta and M_1 are NOT the same as in BWC10 despite the similar names!
	double logeps;	//-1.777 +0.133 -0.146
	double logM_1;	//11.514 +0.053 -0.009
	double epsM_1;	//log10(eps*M_1)
	double alpha;	//-1.412 +.0.020 - 0.105
	double delta13;	//3.508 +0.087 -0.105
	double gamma13;	//0.316 +0.076 -0.012

	double Halomass;

	HaloParams(): beta(0.44), delta(0.57), gamm(1.56), M_1(12.35), M_ss0(5.248e10),
		logeps(-1.777), logM_1(11.514), alpha(-1.412), delta13(3.508), gamma13(0.316), Halomass(0)
	{epsM_1 = log10(pow(10, logeps) * pow(10, logM_1));}
} HP;

double StellarMass(double halomass)
{
	//Return stellar mass of galaxy based on BCW13
	static double S0 = StellarMassFitFunc(0);
	return pow(10, HP.epsM_1 + StellarMassFitFunc(log10(halomass) - HP.logM_1) - S0);
}

double StellarMassFitFunc(double y)
{
	//This equation hurts my brain
	return -log10(pow(10, HP.alpha * y) + 1) + HP.delta13 * pow(log10(1+exp(y)), HP.gamma13)/(1+exp(pow(10, -y)));
}


double HaloMass(double stellarmass)
{
	//Return halo mass of galaxy based on BWC10
	return pow(10, HP.M_1 + HP.beta*log10(stellarmass/HP.M_ss0) + pow(stellarmass/HP.M_ss0, HP.delta)/(1+pow(stellarmass/HP.M_ss0, -HP.gamm)) - 0.5);
}

#endif /* HALOMASS_HPP_ */
