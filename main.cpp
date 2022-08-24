#include "halomass.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>
#include <time.h>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/multi_array.hpp>

using namespace std;
typedef vector<double>::iterator iter;
typedef boost::multi_array<double, 4> FourDArray;

enum GaussIndex {Gamplitude, Goffset, Gsigma};

//Free parameters in the model
struct FreeParameters
{
	FreeParameters(double t, double o, double x, double p): tau0(t), offset(o), x_min(x), texp(p) {}
	double tau0;		//Characteristic merger timescale, default value: 0.57
	double offset;		//Slope in assumed linear offset metallicity - mass ratio relation, default value: 0.3
	double x_min;		//Lower mass ratio limit of integrals below which mergers are supposed to be insignificant, default value: 0.3
	double texp;        //Exponent in merger timescale, default value: -0.3
};

//Parameters for bins symmetrically distributed around 0
struct Binning
{
	Binning(double m, double n): max(m), nBins(n) {Delta = 2*max/nBins;}
	double max;
	double nBins;
	double Delta;
};

struct TabulatedMergerRates
{
	TabulatedMergerRates()
	{
		ifstream datain("merger_rates.txt");
		if (!datain.is_open()) {cout << "Unable to open file merger_rates.txt" << endl;return;}
		string line;
		getline(datain, line);
		int num_massratio = atoi(line.c_str());
		getline(datain, line);
		MassRatioMin = atof(line.c_str());
		getline(datain, line);
		MassRatioDelta = (atof(line.c_str()) - MassRatioMin)/num_massratio;
		getline(datain, line);
		int num_mass = atoi(line.c_str());
		getline(datain, line);
		MassMin = atof(line.c_str());
		getline(datain, line);
		MassDelta = (atof(line.c_str()) - MassMin)/num_mass;

		Rate.resize(num_massratio);
		for (int i = 0; i < num_massratio; i++)
		{
			Rate[i].resize(num_mass);
		}

		getline(datain, line);		//Skip header
		int curMassRatio = 0;
		while (getline(datain, line))
		{
			MassRatio.push_back(atof(line.substr(0, 10).c_str()));
			if (curMassRatio == 0) {Mass.push_back(atof(line.substr(11, 6).c_str()));}
			Rate[curMassRatio][0] =  atof(line.substr(18).c_str());
			for (int curMass = 1; curMass < num_mass; curMass++)
			{
				getline(datain, line);
				if (curMassRatio == 0) {Mass.push_back(atof(line.substr(11, 6).c_str()));}
				Rate[curMassRatio][curMass] =  atof(line.substr(18).c_str());
			}
			curMassRatio++;
		}
	}

	vector<double> MassRatio;
	vector<double> Mass;
	vector<vector <double> > Rate;
	double MassRatioMin;
	double MassRatioDelta;
	double MassMin;
	double MassDelta;
} MergerRate;

struct LambdaCDM
{
	//Parameters for a Lambda-CDM cosmology with zero curvature from Planck 2013 results (http://arxiv.org/pdf/1303.5076v1.pdf)
	//Changed to WMAP 7 cosmology instead to be in accordance with MPA-JHU data
	LambdaCDM(): H_0(70.0), omega_M(0.3), omega_l(0.7) {}
	double H_0;
	double omega_M;
	double omega_l;
} LCDM;

bool GetFMRCoefficients(double (&x)[8], string);
bool GetGaussAndMergerCoefficients(double (&x)[3], double (&y)[5], string);
bool GetOffsetData(string, Binning*, vector<double>&, vector<double>&, string, double x[]);
double ComputeDeviance(vector<double>&, vector<double>&, int);
void DevianceHistogram(vector<double>&, vector<double>&, vector<double>&, int);
void RunModel(FreeParameters*, Binning*, vector<double>&, vector<double>&, double x[], double y[], double);


class Mergers
{
	FreeParameters* Params;
	Binning* ZBins;
	double t_0;
        bool IntegrateTime;

	double TimeToRedshift(double);
	double dMergers(double, double, double);
	double dMergers(double, double);
	double MergerDuration(double);

public:
        Mergers(FreeParameters* P, Binning* B, double t, bool I/*, double Gsigma, double Goffset*/): Params(P), ZBins(B), t_0(t), IntegrateTime(I)
	{}

	double dPdZ(double, double, double x[]);
	double NMergers(double, double, double, double);
	double dMergersTimeIntegrated(double, double,  double);
};

double Mergers::TimeToRedshift(double time)
{
	//Use low-redshift linear approximation for speed
	return 0.071*time+0.0045*time*time;//for Planck cosmology: 0.069*time+0.0006*time*time;
}

//For halo merger rate using halo mass ratio
double Mergers::dMergers(double massratio, double HaloMassTerm, double time)
{
	//Calculate galaxy merger rate according to Fakhouri et al. 2010
	static double InverseHubbleTimeGyr = (LCDM.H_0 * 1000/3.08567758e22)*3.15576e16; //Hubble parameter in units of Gyr^-1
	double z = TimeToRedshift(time);
	double dNdxdz = 0.0104 * HaloMassTerm * pow(massratio, -1.995)*exp(pow(massratio/9.72e-3,0.263));
	//Convert dN/dxdz to dN/dxdt by dN/dxdt=dN/dxdz * dz/dt and dz/dt = H(z)*(1+z) (frequency in Gyr^-1)
	//(1+0.51*z) is an approximation of sqrt(LCDM.omega_M*pow(1+z, 3)+LCDM.omega_l) for omega_M=0.315 and omega_l=0.685
	//(1+0.48*z) is an approximation of sqrt(LCDM.omega_M*pow(1+z, 3)+LCDM.omega_l) for omega_M=0.3 and omega_l=0.7
	//These are valid for z < 0.2
	//If generality is prefered to speed of calculation replace this by sqrt(LCDM.omega_M*pow(1+z, 3)+LCDM.omega_l)
	double dNdxdt = dNdxdz * (1+z) * InverseHubbleTimeGyr * (1+0.48*z);//(1+0.51*z);
	return dNdxdt;
}

//For galaxy merger rate using stellar mass ratio
//From lookup table using bilinear interpolation between given mass ratios and stellar masses
//Based on halo merger rate of Fakhouri et al. 2010 and the stellar mass-halo mass relation of Behroozi et al. 2010
//Assume z=0.1 and a H_0=70, omega_M=0.3, omega_l=0.7 cosmology
double Mergers::dMergers(double StellarMassRatio, double StellarMass)
{
	int x1pos = int((StellarMassRatio-MergerRate.MassRatioMin)/MergerRate.MassRatioDelta);
	if (MergerRate.MassRatio[x1pos] > StellarMassRatio) {x1pos--;}
	double x1 = MergerRate.MassRatio[x1pos];
	double y1 = MergerRate.Mass[int((StellarMass-MergerRate.MassMin)/MergerRate.MassDelta)];
	double x2 = MergerRate.MassRatio[x1pos + 1];
	double y2 = MergerRate.Mass[int((StellarMass-MergerRate.MassMin)/MergerRate.MassDelta) + 1];
	double f11 = MergerRate.Rate[x1pos][int((StellarMass-MergerRate.MassMin)/MergerRate.MassDelta)];
	double f12 = MergerRate.Rate[x1pos][int((StellarMass-MergerRate.MassMin)/MergerRate.MassDelta) + 1];
	double f21 = MergerRate.Rate[x1pos + 1][int((StellarMass-MergerRate.MassMin)/MergerRate.MassDelta)];
	double f22 = MergerRate.Rate[x1pos + 1][int((StellarMass-MergerRate.MassMin)/MergerRate.MassDelta) + 1];

	double x = StellarMassRatio;
	double y = StellarMass;
	double dNdxdt = (1/((x2 - x1)*(y2 - y1))) * (f11*(x2 - x)*(y2 - y) + f21*(x - x1)*(y2 - y) + f12*(x2 - x)*(y - y1) + f22*(x - x1)*(y - y1));
	return dNdxdt;
}

double Mergers::dPdZ(double massratio, double Z, double GaussCoeffs[])
{
	static double sqrt2pi = 2.50662827;

	double A = 1/(GaussCoeffs[Gsigma]*sqrt2pi);
	return A*exp(-(Z - (GaussCoeffs[Goffset] - Params->offset * massratio))*(Z - (GaussCoeffs[Goffset] - Params->offset * massratio))/(2*GaussCoeffs[Gsigma]*GaussCoeffs[Gsigma]));
}

double Mergers::MergerDuration(double massratio)
{
	//Time that the merger is non-negligible in Gyr
	//using n2med, first passage timescale from http://adsabs.harvard.edu/abs/2008MNRAS.384..386C
	if (Params->texp == 0.0) {return Params->tau0;}
	else if (Params->texp == -1.0) {return Params->tau0 / massratio;}
	else {return Params->tau0 * pow(massratio, Params->texp);}
	//return Params->tau0;
}

double Mergers::dMergersTimeIntegrated(double Delta_t, double Mass, double x)
{
	double t_span = MergerDuration(x);		//Span of time in Gyr
	double dNmergerdx = 0;
	dNmergerdx = dMergers(x, Mass) *  t_span;
	return dNmergerdx;
}

double Mergers::NMergers(double Delta_t, double x_min, double Delta_x, double Mass)
{
	double R = 0;
	for(double x = x_min; x + Delta_x <= 1; x += Delta_x)
	{
                double R1 = dMergersTimeIntegrated(Delta_t, Mass, x);
                double R2 = dMergersTimeIntegrated(Delta_t, Mass, x + Delta_x);
		R += 0.5 * (R1 + R2) * Delta_x;
	}
	return R;
}

int main(int argc, char* argv[])
{
	string Z_maxstr = argv[1];
	string nBinsstr = argv[2];
	string tau_minstr = argv[3];
	string tau_maxstr = argv[4];
	string delta_taustr = argv[5];
	string offset_minstr = argv[6];
	string offset_maxstr = argv[7];
	string delta_offsetstr = argv[8];
	string x_min_minstr = argv[9];
	string x_min_maxstr = argv[10];
	string delta_x_minstr = argv[11];
    string texp_minstr = argv[12];
    string texp_maxstr = argv[13];
    string delta_texpstr = argv[14];
	string datafile = argv[15];
	string FMRtype = argv[16];

	const double Z_max = atof(Z_maxstr.c_str());
	int nBins = atoi(nBinsstr.c_str());
	if ((nBins < 1) || (Z_max < 0))
	 {cout << "Maximum Z must be non-negative and the number of bins must be at least 1!";return -1;}
	Binning ZBins(Z_max, nBins);
	cout << nBins << " metallicity offset bins from " << -Z_max << " to " << Z_max << " spaced by " << ZBins.Delta << endl;

	const double tau0_min = atof(tau_minstr.c_str());
	const double tau0_max = atof(tau_maxstr.c_str());
	const double delta_tau0 = atof(delta_taustr.c_str());
	const double offset_min = atof(offset_minstr.c_str());
	const double offset_max = atof(offset_maxstr.c_str());
	const double delta_offset = atof(delta_offsetstr.c_str());
	const double x_min_min = atof(x_min_minstr.c_str());
	const double x_min_max = atof(x_min_maxstr.c_str());
	const double delta_x_min = atof(delta_x_minstr.c_str());
    const double texp_min = atof(texp_minstr.c_str());
    const double texp_max = atof(texp_maxstr.c_str());
    const double delta_texp = atof(delta_texpstr.c_str());

	int num_tau0 = (tau0_max - tau0_min)/delta_tau0 + 1;
	int num_offset = (offset_max - offset_min)/delta_offset + 1;
	int num_x_min = (x_min_max - x_min_min)/delta_x_min + 1;
	int num_texp = (texp_max - texp_min)/delta_texp + 1;

	vector<double> MetalOffset(nBins);
	vector<double> MetalOffsetDist_model(ZBins.nBins, 0);
	vector<double> MetalOffsetDist_data(ZBins.nBins, 0);

	double GaussCoeffs[3];
	double MergerWeights[5];
	if (!GetGaussAndMergerCoefficients(GaussCoeffs, MergerWeights, FMRtype)) {cout << "Failed to load coeffiecients for Gaussian - Quitting." << endl;return -1;}

	FourDArray::extent_gen extents;
	FourDArray Deviance(extents[num_tau0][num_offset][num_x_min][num_texp]);	//Deviance for models with parameters (tau0, offset, x_min, texp) compared to data

	for (int BinZ = 0; BinZ < nBins; BinZ++)
	  {
	    MetalOffset[BinZ] = -ZBins.max + BinZ*ZBins.Delta;
	  }

	if (!GetOffsetData(datafile, &ZBins, MetalOffset, MetalOffsetDist_data, FMRtype, GaussCoeffs)) {cout << "Failed to load metallicity data - Quitting." << endl;return -1;}

	int DevianceUpperLimBin = nBins - 1;
	double SumWithHighZTail = 0;
	double SumExcludingHighZTail = 0;
	bool LimFound = false;

	for(int BinZ = 0; BinZ < nBins; BinZ++)
	{
		if ((MetalOffset[BinZ] > GaussCoeffs[Goffset] + 2*GaussCoeffs[Gsigma]) && (!LimFound))
		{
			DevianceUpperLimBin = BinZ;
			LimFound = true;
		}
		if (!LimFound) {SumExcludingHighZTail += MetalOffsetDist_data[BinZ];}
		SumWithHighZTail += MetalOffsetDist_data[BinZ];
	}

	//Coefficient to correct for high Z tail not being included in model when renormalizing model from probability to counts
	double HighZTailCorrection = SumExcludingHighZTail/SumWithHighZTail;

	cout << num_tau0 << " tau0 values, " << num_offset << " offset values, " << num_x_min << " x_min values, " << num_texp << " texp values" << endl;
	cout << "Creating models and calculating log-likelihood values:" << endl;

	//Run outer loop on multiple threads if possible. Note that this cannot take advantage of more than num_tau0 threads.
	#pragma omp parallel for firstprivate(MetalOffsetDist_model, MetalOffsetDist_data)
	for (int i = 0; i < num_tau0; ++i)
	{
		double tau0 = tau0_min + i*delta_tau0;
		#ifdef _OPENMP
			if (omp_get_num_threads() == 1) {cout << (i/float(num_tau0))*100 << "%" << endl;}
		#else
			cout << (i/float(num_tau0))*100 << "%" << endl;
		#endif
		for (int j = 0; j < num_offset; ++j)
		{
			double offset = offset_min + j*delta_offset;
			for (int k = 0; k < num_x_min; ++k)
			{
				double x_min = x_min_min + k*delta_x_min;
				for (int l = 0; l < num_texp; ++l)
				{
					double texp = texp_min + l*delta_texp;
					fill(MetalOffsetDist_model.begin(), MetalOffsetDist_model.end(), 0.0);
					FreeParameters CurParams(tau0, offset, x_min, texp);
					RunModel(&CurParams, &ZBins, MetalOffset, MetalOffsetDist_model, GaussCoeffs, MergerWeights, HighZTailCorrection);
                                        Deviance[i][j][k][l] = ComputeDeviance(MetalOffsetDist_data, MetalOffsetDist_model, DevianceUpperLimBin);
				}
			}
		}
		#ifdef _OPENMP
		if (omp_get_num_threads() > 1) {cout << "tau0 = " << tau0 << " models calculated by thread " << omp_get_thread_num() << endl;}
		#endif
	}
	cout << "100%" << endl;

	string filename = "likelihood.txt";
	ofstream dataout(filename.c_str());
	if (!dataout.is_open()) {cout << "Unable to open file likelihood.txt for output" << endl;return -1;}
	dataout << "Arguments:\t";
	for (int argn = 1; argn < 16; argn++)
	{
		dataout << argv[argn] << " ";
	}
	dataout << argv[16] << endl << "tau0\toffset\tx_min\ttexp\t-log likelihood (with constant term sum(ln(d_i!)) dropped)" << endl;

	//Find smallest negative log-likelihood and best-fit parameters
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	double DevianceMin = 1e70;
	vector<FreeParameters> DevianceMinPos;
	FreeParameters Blank(-1,-1,-1,-1);
	DevianceMinPos.push_back(Blank);
	for (double tau0 = tau0_min; tau0 <= tau0_max; tau0 += delta_tau0)
	{
		j = 0;
		for (double offset = offset_min; offset < offset_max; offset += delta_offset)
		{
			k = 0;
			for (double x_min = x_min_min; x_min < x_min_max; x_min += delta_x_min)
			{
				l = 0;
				for (double texp = texp_min; texp < texp_max; texp += delta_texp)
				{
					dataout << setprecision(10) << tau0 << "\t" << offset << "\t" << x_min << "\t" << texp << "\t" << Deviance[i][j][k][l] << endl;
					if (Deviance[i][j][k][l] <= DevianceMin)
					{
						if (Deviance[i][j][k][l] < DevianceMin) {DevianceMinPos.clear();DevianceMin = Deviance[i][j][k][l];}
						DevianceMinPos.push_back(FreeParameters(tau0, offset, x_min, texp));
					}
					l++;
				}
				k++;
			}
			j++;
		}
		i++;
	}

        cout << "---------------------------------------------" << endl;
        cout << "Minimum shifted log-likelihood:\t" << DevianceMin << endl << "Minimum log-likelihood parameters:" << endl;

	dataout << endl << "Minimum shifted log-likelihood:\t" << DevianceMin << endl << "Minimum log-likelihood parameters:" << endl;
	for (unsigned int i = 0; i < DevianceMinPos.size(); i++)
	{
		dataout << "Min logL model no. " << i+1 << ": tau0=" << DevianceMinPos[i].tau0 << " offset=" << DevianceMinPos[i].offset << " x_min=" << DevianceMinPos[i].x_min << " texp=" << DevianceMinPos[i].texp << endl;
		cout << "tau0=" << DevianceMinPos[i].tau0 << " offset=" << DevianceMinPos[i].offset << " x_min=" << DevianceMinPos[i].x_min << " texp=" << DevianceMinPos[i].texp << endl;
	}
	dataout.close();

	cout << "---------------------------------------------" << endl;
	cout << "Writing best-fit model metallicity distribution to file" << endl;
	fill(MetalOffsetDist_model.begin(), MetalOffsetDist_model.end(), 0.0);
	FreeParameters CurParams = DevianceMinPos[0];
	RunModel(&CurParams, &ZBins, MetalOffset, MetalOffsetDist_model, GaussCoeffs, MergerWeights, HighZTailCorrection);
	vector<double> DevianceHist;
	DevianceHistogram(DevianceHist, MetalOffsetDist_data, MetalOffsetDist_model, DevianceUpperLimBin);

	filename = "deviance_histogram.txt";
	ofstream histout(filename.c_str());
	if (!histout.is_open()) {cout << "Unable to open file deviance_histogram.txt for output" << endl;return -1;}
	histout << "Zoffset\tDeviance" << endl;
	for (unsigned int i = 0; i < DevianceHist.size(); i++)
	{
		histout << MetalOffset[i] << "\t" << DevianceHist[i] << endl;
	}
	histout.close();

	filename = "MetalOffsetDist_model.txt";
	ofstream metalout(filename.c_str());
	if (!metalout.is_open()) {cout << "Unable to open file MetalOffsetDist_model.txt for output" << endl;return -1;}
	metalout << "Zoffset\tProbability" << endl;
	for (int i = 0; i < nBins; i++)
	{
		metalout << MetalOffset[i] << "\t" << MetalOffsetDist_model[i] << endl;
	}
	metalout.close();
	cout << "---------------------------------------------" << endl;

	return 0;
}

//Read FMR/MZR coefficients from file. Form: c[0] + c[1]*m + c[2]*m^2 + c[3]*m^3 + c[4]*m^4 + c[5]*s + c[6]*s^2 + c[7]*m*s
bool GetFMRCoefficients(double (&Coeffs)[8], string Type)
{
	ifstream datain("FMRCoeffs.txt");
	if (!datain.is_open()) {cout << "Unable to open file FMRCoeffs.txt"<< endl;return false;}
	string line;
	int i = 0;
	bool typefound = false;
	while (getline(datain, line))
	{
		if (line == "Type: " + Type) {typefound=true;continue;}
		if (!typefound) {continue;}
		Coeffs[i] = atof(line.c_str());
		if (i == 7) {break;}
		i++;
	}
	datain.close();
	if (i != 7) {return false;} else {return true;}
}

//Read Gaussian coefficients and merger mass weights from file. Form: c[0]=amplitude, c[1]=mean, c[2]=scatter
bool GetGaussAndMergerCoefficients(double (&Coeffs)[3], double (&Weights)[5], string Type)
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
			if (i < 3) {Coeffs[i] = atof(line.c_str());}
			else {Weights[i-3] = atof(line.c_str());}
			if (i == 2) {getline(datain, line);i++;continue;}
			if (i == 7) {break;}
			i++;
		}
	}
	datain.close();
	if (i != 7) {return false;} else {return true;}
}

//Read metallicity data from file and calculate offset from FMR
bool GetOffsetData(string filename, Binning* ZBins, vector<double>& MetalOffset, vector<double>& MetalOffsetDist_data, string Type, double GaussCoeffs[])
{
	//Structure of data file: GalaxyID log(O/H)+12 log(M_s) log(SFR) Redshift, header at first line
	double FMRCoeffs[8];		//Coeffecients for FMR/MZR
	if (!GetFMRCoefficients(FMRCoeffs, Type)) {return false;}
	ifstream datain(filename.c_str());
	if (!datain.is_open()) {cout << "Unable to open file " << filename << endl;return false;}
	vector<double> FMROffset;
	string line;
	getline(datain, line);		//Skip header
	while (getline(datain, line))
	{
		int OHbegin = line.find(" ")+1;
		int OHend = line.find(" ",OHbegin+1)-1;
		int Mend = line.find(" ", OHend+2)-1;
		int SFRend = line.find(" ", Mend+2)-1;
		double OH = atof(line.substr(OHbegin,OHend - OHbegin + 1).c_str());		//Second column of data file: log(O/H)+12
		double m = atof(line.substr(OHend+1, Mend - OHend + 1).c_str()) - 10;	//Third column of data file: log(Stellar mass), subtract 10 to get m parameter in FMR
		double s = atof(line.substr(Mend+2, SFRend - Mend + 1).c_str());		//Fourth column of data file: logSFR, s parameter in FMR

		FMROffset.push_back(OH - (FMRCoeffs[0] + FMRCoeffs[1]*m + FMRCoeffs[2]*m*m + FMRCoeffs[3]*m*m*m + FMRCoeffs[4]*m*m*m*m + FMRCoeffs[5]*s + FMRCoeffs[6]*s*s + FMRCoeffs[7]*m*s));
	}
	int OutofBounds = 0;
	int sum = 0;
	for(unsigned int i = 0; i < FMROffset.size(); i++)
	{
		if ((FMROffset[i] < -ZBins->max) || (FMROffset[i] > ZBins->max))
                {
			OutofBounds++;
		}
		else
		{
			for (int curBin = 1; curBin < ZBins->nBins; curBin++)
			{
			  if (FMROffset[i] < MetalOffset[curBin]) {MetalOffsetDist_data[curBin-1]++;sum++;break;}
			}
		}
	}
	cout << "Out of Bounds: " << OutofBounds << " of " << sum << endl;

	ofstream dataout("MetalOffsetDist_data.txt");
	if (!dataout.is_open()) {cout << "Unable to open file MetalOffsetDist_data.txt for output" << endl;return -1;}
	dataout << "Zoffset\tProbability" << endl;
	for (int i = 0; i < ZBins->nBins; i++)
	{
		MetalOffset[i] = MetalOffset[i] + 0.5*ZBins->Delta;
		dataout << MetalOffset[i] << "\t" << MetalOffsetDist_data[i] << endl;
	}
	dataout.close();

	return true;
}

//Compute -log likelihood with the constant ln(d!) term dropped, assuming Poisson distributed errors in each bin
double ComputeDeviance(vector<double>& MetalOffsetDist_data, vector<double>& MetalOffsetDist_model, int DevianceUpperLimBin)
{
	double Dev = 0;
	for (int m = 0; m <= DevianceUpperLimBin; m++)
	{
		Dev += MetalOffsetDist_model[m] - MetalOffsetDist_data[m]*log(MetalOffsetDist_model[m]);
	}
	return Dev;//2*Dev;
}

//Compute the deviance contribution for each bin
//The deviance is -2*log likelihood ratio of model given data to data given data assuming Poisson distributed errors in each bin
void DevianceHistogram(vector<double>& DevianceHist, vector<double>& MetalOffsetDist_data, vector<double>& MetalOffsetDist_model, int DevianceUpperLimBin)
{
	for (int m = 0; m <= DevianceUpperLimBin; m++)
	{
		if (MetalOffsetDist_data[m] != 0)
		{
			DevianceHist.push_back(2*(MetalOffsetDist_data[m]*log(MetalOffsetDist_data[m]/MetalOffsetDist_model[m]) - (MetalOffsetDist_data[m] - MetalOffsetDist_model[m])));
		}
		else
		{
			DevianceHist.push_back(2*(MetalOffsetDist_data[m] - MetalOffsetDist_model[m]));
		}
	}
}


void RunModel(FreeParameters* Params, Binning* ZBins, vector<double>& MetalOffset, vector<double>& MetalOffsetDist_model, double GaussCoeffs[], double MergerWeights[], double HighZTailCorrection)
{
	//const double t_now = 13.811;
	double Delta_t = 0.001;
	//double Delta_x = 0.01;
	const double Delta_x_fine = 1e-4;
	Mergers Merger(Params, ZBins, 0, false);
	double Nnew = 0;

	int i = 0;
	for (double M_sb = 9.25;M_sb < 11.5; M_sb += 0.5)
	{
		Nnew += Merger.NMergers(Delta_t, Params->x_min, Delta_x_fine, M_sb/*HaloMassTerm*/) * MergerWeights[i];

		double x_init = 0;
		x_init=Params->x_min;
		for (double x = x_init; x <= 1; x += Delta_x_fine) //loop over mass ratio bins
		{
                        double dNmergerdx = Merger.dMergersTimeIntegrated(Delta_t, M_sb, x);
			for(int BinZ = 0; BinZ < ZBins->nBins; BinZ++)
			{
				double curZ = -ZBins->max + BinZ*ZBins->Delta;
				MetalOffsetDist_model[BinZ] += (dNmergerdx * Merger.dPdZ(x, curZ, GaussCoeffs) * Delta_x_fine) * MergerWeights[i];
			}
		}
		i++;
	}

	//Assume that the probability distribution of DeltaZ in absence of mergers is a Gaussian and add this to the merger metallicity distribution
	double sum = 0;
	for (int i = 0; i < ZBins->nBins; i++)
	{
		MetalOffsetDist_model[i] += (1 - Nnew) * (1./(sqrt(2*M_PI)*GaussCoeffs[Gsigma])) * exp(-((MetalOffset[i]-GaussCoeffs[Goffset])*(MetalOffset[i]-GaussCoeffs[Goffset]))/(2*GaussCoeffs[Gsigma]*GaussCoeffs[Gsigma]));
	    sum += MetalOffsetDist_model[i] * ZBins->Delta * HighZTailCorrection;
	}
	double sum_after = 0;

	//Renormalize distribution from probability to galaxy counts
	for (iter j = MetalOffsetDist_model.begin();j != MetalOffsetDist_model.end();++j)
	{
		*j *= GaussCoeffs[Gsigma]*sqrt(2*M_PI)*GaussCoeffs[Gamplitude]/sum;
		sum_after += *j;
	}
}
