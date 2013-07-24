/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#include "interface/Constants.h"

namespace Elegent
{

//----------------------------------------------------------------------------------------------------

TComplex i(0, 1);
Constants *cnts = NULL;

// physics constants
double Constants::alpha =			7.297E-3;
double Constants::proton_mass =		0.938271;
double Constants::neutron_mass =	0.939565;
double Constants::hbarc =			0.197326;
double Constants::sq_hbarc =		0.389379;
double Constants::M = proton_mass;
double Constants::M_sq= proton_mass * proton_mass;
double Constants::kappa =			1.793;

// mathematical constants
double Constants::pi =				M_PI;
double Constants::gamma =			0.577215;

//----------------------------------------------------------------------------------------------------

Constants::Constants(double W, Constants::ParticleMode mode)
{
	sqrt_s = W;											/// sqrt_s / GeV
	s = sqrt_s*sqrt_s;									/// s / GeV^2
	ln_s = log(s);										/// ln(s / GeV^2)
	p_cms = sqrt(s /4. - proton_mass * proton_mass);	/// particle CMS impuls
	sig_fac = sq_hbarc * pi / (s * p_cms * p_cms);		/// d sig/dt = sig_fac * |A|^2
	t_min = 4.*proton_mass * proton_mass - s;			/// NEGATIVE !
	
	pMode = mode;
	//printf("Running at \e[00;34mW = %.1f GeV, in ", sqrt_s);
	//if (pMode == mPP) printf("p-p");
	//if (pMode == mAPP) printf("p-anti p");
	//printf(" mode\e[00m\n");
}

//----------------------------------------------------------------------------------------------------

void Constants::Init(double W, Constants::ParticleMode mode)
{
	cnts = new Constants(W, mode);
}

//----------------------------------------------------------------------------------------------------

void Constants::Print()
{
	printf(">> Constants::Print\n");
	printf("	alpha = %E\n", alpha);
	printf("	proton_mass = %E\n", proton_mass);
	printf("	neutron_mass = %E\n", neutron_mass);
	printf("	hbarc = %E\n", hbarc);
	printf("	sq_hbarc = %E\n", sq_hbarc);
	printf("	M = %E\n", M);
	printf("	M_sq = %E\n", M_sq);
	printf("	pi = %E\n", pi);
	printf("	gamma = %E\n", gamma);
	printf("	sqrt_s = %E\n", sqrt_s);
	printf("	s = %E\n", s);
	printf("	ln_s = %E\n", ln_s);
	printf("	p_cms = %E\n", p_cms);
	printf("	sig_fac = %E\n", sig_fac);
	printf("	t_min = %E\n", t_min);
	printf("	pMode = %s\n", (pMode == mPP) ? "pp" : "app");
}

} // namespace
