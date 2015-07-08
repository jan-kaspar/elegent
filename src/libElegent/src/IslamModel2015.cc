/********************************************************************************

    Copyright 2013 Jan Kašpar

    This file is part of Elegent (http://elegent.hepforge.org/).

    Elegent is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Elegent is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Elegent.  If not, see <http://www.gnu.org/licenses/>.
 
********************************************************************************/

#include "interface/IslamModel2015.h"

#include "interface/Constants.h"

using namespace std;
using namespace Elegent;


//----------------------------------------------------------------------------------------------------

IslamModel2015::IslamModel2015()
{
	fullLabel.name = "Islam et al."; shortLabel.name = "islam";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

IslamModel2015::~IslamModel2015()
{
	if (integ_workspace_initialized)
	{
		gsl_integration_workspace_free(integ_workspace_b);
		gsl_integration_workspace_free(integ_workspace_t);
	}
}

//----------------------------------------------------------------------------------------------------

void IslamModel2015::Configure(IslamModel2015::ModeType _m)
{
	mode = _m;
	
	fullLabel.variant = ""; shortLabel.variant = "";
	fullLabel.version = "EDS'15"; shortLabel.version = "15";

	if (mode == mDiff)
		{ fullLabel.mode = "diffraction"; shortLabel.mode = "diff"; }
	if (mode == mOmega)
		{ fullLabel.mode = "omega"; shortLabel.mode = "omega"; }
	if (mode == mQuark)
		{ fullLabel.mode = "quark-quark"; shortLabel.mode = "quark"; }
	if (mode == mPolarization)
		{ fullLabel.mode = "polarization"; shortLabel.mode = "pol"; }
	if (mode == mFull)
		{ fullLabel.mode = "full"; shortLabel.mode = "full"; }
}

//----------------------------------------------------------------------------------------------------

void IslamModel2015::Init()
{
	// ---------- diffraction amplitude ----------
	// [2]
	R0 = 3.23;
   	R1 = 0.0508;
	a0 = 0.375;
	a1 = 0.109;
	R = TComplex(R0 + R1*cnts->ln_s, -R1*cnts->pi/2.);
	a = TComplex(a0 + a1*cnts->ln_s, -a1*cnts->pi/2.);
	
	// [2]
	eta0 = 0.0594;
	c0 = 0.0;
	sigma = 1.0;
	
	// function g(s) according to Eq. (4.2) in [2006 review paper]
	g_s = (1. - eta0 + c0 / TComplex::Power(-i * cnts->s, sigma)) * (1. + TComplex::Exp(-R/a)) / (1. - TComplex::Exp(-R/a));
	
	// ---------- screening factor due to diffraction ----------
	
	// [2]
	lambda0 = 1.23;
	d0 = 10.6;
	s0 = 1000. * 1000.;	// GeV^2

	double sign = (cnts->pMode == cnts->mPP) ? +1. : -1.;
	fac_scr_diff = (eta0 + c0 / TComplex::Power(-i * cnts->s, sigma))
		+ sign * i * (lambda0 - d0 / pow(cnts->s/s0, 2.));
	
	// ---------- omega-exchange amplitude ----------

	// [2]
	gamma_hat = 1.94;
	theta_hat = 1.366;
	beta = 3.075;
	m_omega = 0.801;

	// ---------- screening due to omega exchange ----------

	// [2]
	b_tilde = 0.;

	fac_scr_omega = TComplex::Exp(i * Chi_omega(b_tilde));
	
	// ---------- quark-quark amplitude ----------

	// [2]
	gamma_gg = 0.101;
	lambda = 0.221;
	m = 1.67;
	mu = 1./4.;

	m0 = 3.464;
	M = 0.939;
	
	// ---------- polarisation amplitude ----------

	ht = 0.024;
	wd = 0.0034;
	sf = 0.647;
	
	// integration parameters
	upper_bound_b = 50.;
	precision_b = 5E-2;

	upper_bound_t = -15.;
   	precision_t = 1E-3;

	// prepare integration workspace
	if (!integ_workspace_initialized)
	{
		integ_workspace_size_b = 100;
		integ_workspace_b = gsl_integration_workspace_alloc(integ_workspace_size_b);

		integ_workspace_size_t = 100;
		integ_workspace_t = gsl_integration_workspace_alloc(integ_workspace_size_b);

		integ_workspace_initialized = true;
	}
}

//----------------------------------------------------------------------------------------------------

void IslamModel2015::Print() const
{
	printf(">> IslamModel2015::Print\n");
	printf("\t%s\n", CompileFullLabel().c_str());

	printf("\tdiffraction amplitude\n");
	printf("\t\tR0=%f, R1=%f\n", R0, R1);
	printf("\t\ta0=%f, a1=%f\n", a0, a1);
	printf("\t\teta0=%f, c0=%f, sigma=%f\n", eta0, c0, sigma);

	printf("\tscreening due to diffraction\n");
	printf("\t\tlambda0=%f, d0=%f, s0=%E\n", lambda0, d0, s0);

	printf("\tomega-exchange amplitude\n");
	printf("\t\tgamma_hat=%f, theta_hat=%f, beta=%f, m_omega=%f\n", gamma_hat, theta_hat, beta, m_omega);

	printf("\tscreening due to omega exchange\n");
	printf("\t\tb_tilde=%f\n", b_tilde);

	printf("\tquark confinement\n");
	printf("\t\tm0=%f, M=%f\n", m0, M);

	printf("\tquark-quard scattering variables\n");
	printf("\t\tgamma_gg=%f, lambda=%f, m=%f, mu=%f\n", gamma_gg, lambda, m, mu);

	printf("\tpolarisation amplitude\n");
	printf("\t\tht=%f, wd=%f, sf=%f\n", ht, wd, sf);

	printf("\tintegration variables\n");
	printf("\t\tprecision_b = %E\n", precision_b);
	printf("\t\tprecision_t = %E\n", precision_t);
	printf("\t\tupper_bound_b = %E\n", upper_bound_b);
	printf("\t\tupper_bound_t = %E\n", upper_bound_t);
}

//-------------------------------------- DIFFRACTION AMPLITUDE ------------------------------------

TComplex IslamModel2015::GammaDPlus(double b) const
{
	/// Eq. (2) in [1], without g(s)
	return 1. / (1. + TComplex::Exp((b - R) / a)) + 1. / (1. + TComplex::Exp((-b - R) / a)) - 1.;
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel2015::GammaDPlus_J0(double b, double *par, const void *vobj)
{
	const IslamModel2015 *obj = (IslamModel2015 *) vobj;
	const double &t = par[0];	

	return obj->GammaDPlus(b) * b * TMath::BesselJ0(b*sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel2015::T_diff(double t) const
{
	/// t < 0
	double par[] = { t };
	TComplex scale = g_s * i * cnts->p_cms * cnts->sqrt_s;
	return scale * ComplexIntegrate(GammaDPlus_J0, par, this, 0., upper_bound_b, 0., precision_b,
		integ_workspace_size_b, integ_workspace_b, "IslamModel2015::T_diff");
}


//----------------------------------------- CORE AMPLITUDE ----------------------------------------

TComplex IslamModel2015::Chi_omega(double b) const
{
	double sign = (cnts->pMode == cnts->mPP) ? +1. : -1.;

	return -sign * 2. * gamma_hat * TComplex::Exp(sign * i * theta_hat) * TMath::BesselK0(m_omega * sqrt(b*b + beta*beta));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel2015::GammaOmega_J0(double b, double *par, const void *vobj)
{
	const IslamModel2015 *obj = (IslamModel2015 *) vobj;
	const double &t = par[0];	
	
	return b * TMath::BesselJ0(b*sqrt(-t)) * (1. - TComplex::Exp(i * obj->Chi_omega(b)));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel2015::T_omega(double t) const
{
	double par[] = { t };
	TComplex scale = fac_scr_diff * i * cnts->p_cms * cnts->sqrt_s;
	return scale * ComplexIntegrate(GammaOmega_J0, par, this, 0., upper_bound_b, 0., precision_b,
		integ_workspace_size_b, integ_workspace_b, "IslamModel2015::T_omega");
}

//---------------------------------------- QUARK AMPLITUDE ----------------------------------------

double IslamModel2015::I_integral(double qt, double al) const
{
	double ap = qt/2./al;
	double a = sqrt(ap * ap + 1.);
	double ival;
	
	// for small qt values just substitute limit value 16/3
	if (qt > 1E-10)
		ival = (	2./a/a + 1./ap/ap - 3.*ap*ap/a/a/a/a	) / a/ap * log(a + ap)	-	1./a/a/ap/ap + 3./a/a/a/a;
	else ival = 16./3.;

	return 1./8. /al/al/al/al * ival;
}

//----------------------------------------------------------------------------------------------------

double IslamModel2015::F_cal_integ(double x, double *par, const void *vobj)
{
	IslamModel2015 *obj = (IslamModel2015 *) vobj;

	const double &qt = par[0];
	const double &lambda = par[1];
	const double &m0sq = par[2];

	double al_sq = m0sq/4. + obj->M*obj->M * x * x;
	double al = sqrt(al_sq);

	return pow(x, 1. + lambda) / al_sq * obj->I_integral(qt, al);
}

//----------------------------------------------------------------------------------------------------

double IslamModel2015::F_cal(double qt) const
{
	double par[] = { qt, lambda, m0*m0 };
	double I = RealIntegrate(F_cal_integ, par, this, 0., 1., 0., 1E-3, integ_workspace_size_b,
		integ_workspace_b, "IslamModel2015::F_cal");

	return M * pow(m0, 5.) / 8. / cnts->pi * I;
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel2015::T_quark(double t) const
{
	double q = sqrt(-t);
	double qt = sqrt(-t * (-t / cnts->t_min + 1.));
	double F = F_cal(qt);
	TComplex T_gg = i * cnts->s * gamma_gg * TComplex::Power(-i * cnts->s, lambda) * F*F / pow(1. + q*q/m/m, 2.*(mu+1.));
	return fac_scr_diff * fac_scr_omega * T_gg;
}

//------------------------------------ POLARIZIZATION AMPLITUDE --------------------------------------

TComplex IslamModel2015::T_polarization(double t) const
{
	double sign = (cnts->pMode == cnts->mPP) ? -1. : +1.;
	double q = sqrt(-t);
	return sign * cnts->p_cms * cnts->sqrt_s / 2. / wd * ht * exp(- (sf*sf + q*q) / 4. / wd)
		* TMath::BesselI0(sf * q / 2. / wd);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------- FULL AMPLITUDE -------------------------------------------

TComplex IslamModel2015::Amp(double t) const
{	
	switch (mode)
	{
		case mDiff:			return T_diff(t);
		case mOmega:		return T_omega(t);
		case mQuark:		return T_quark(t);
		case mPolarization:	return T_polarization(t);
		case mFull:	 		return T_diff(t) + T_omega(t) + T_quark(t) + T_polarization(t);

		default:
			printf("ERROR in IslamModel2015::Amp > unknown mode %i\n", mode);
			abort();
	}
}

//----------------------------------------------------------------------------------------------------
//---------------------------------------- PROFILE FUNCTIONS --------------------------------------

TComplex IslamModel2015::Amp_J0(double t, double *par, const void *vobj)
{
	const IslamModel2015 *obj = (IslamModel2015 *) vobj;
	const double &b = par[0];	

	return obj->Amp(t) * TMath::BesselJ0(b*sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel2015::Prf(double b_fm) const
{
	double b = b_fm / cnts->hbarc;	// b in GeV^-1
	double par[] = { b };

	TComplex I = ComplexIntegrate(Amp_J0, par, this, upper_bound_t, 0., 0., precision_t,
		integ_workspace_size_t, integ_workspace_t, "IslamModel2015::Prf");
	return I / 4. / cnts->p_cms / cnts->sqrt_s;
}
