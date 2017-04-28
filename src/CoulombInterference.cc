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

#include "interface/CoulombInterference.h"
#include "interface/Constants.h"
#include "interface/IFunctionInterpolator.h"

using namespace std;
using namespace Elegent;

namespace Elegent
{

//----------------------------------------------------------------------------------------------------

// a global instance
CoulombInterference *coulomb = new CoulombInterference();

//----------------------------------------------------------------------------------------------------

CoulombInterference::CoulombInterference() : mode(mPC), ffType(ffPuckett),
	tau(1E-10), T(10.), precision(1E-5)
{
	integ_workspace_size = 1000;
	integ_workspace = gsl_integration_workspace_alloc(integ_workspace_size);
	integ_workspace2 = gsl_integration_workspace_alloc(integ_workspace_size);
}

//----------------------------------------------------------------------------------------------------

CoulombInterference::~CoulombInterference()
{
	gsl_integration_workspace_free(integ_workspace);
	gsl_integration_workspace_free(integ_workspace2);

	if (iFunctionInterpolator)
		delete iFunctionInterpolator;
}

//----------------------------------------------------------------------------------------------------

void CoulombInterference::Print() const
{
	printf(">> CoulombInterference::Print\n");
	printf("\tmode: %s\n", GetModeString().c_str());
	printf("\tform factor: %s\n", GetFFName().c_str());
	printf("\tT = %E\n", T);
	printf("\ttau = %E\n", tau);
	printf("\tprecision = %E\n", precision);
	printf("\tuseIFunctionInterpolator = %u\n", useIFunctionInterpolator);
}

//----------------------------------------------------------------------------------------------------

string CoulombInterference::GetModeString() const
{
	switch (mode)
	{
		case mPC: return "PC";
		case mPH: return "PH";
		case mWY: return "WY";
		case mSWY: return "SWY";
		case mKL: return "KL";
		case mCahn: return "Cahn";
		default: return "unknown";
	}
}

//--------------------------------------- FORM FACTORS ---------------------------------------------

double di_la_sq = 0.71;

double bor_E_msq[] = {0.13745,	0.58485,	1.7164,	6.0042};
double bor_E_c[] =	 {0.0301,	0.8018,		-1.0882, 0.2642};

double bor_M_msq[] = {0.339, 0.583, 1.711, 13.793};
double bor_M_c[] =	 {0.269, 0.346, -0.672, 0.069};

double kel_E_a[] = {1., -0.24};
double kel_E_b[] = {1., 10.98, 12.82, 21.97};

double kel_M_a[] = {1., 0.12};
double kel_M_b[] = {1., 10.97, 18.86, 6.55};

double arr_E_a[] = {1., 3.439, -1.602, 0.068};
double arr_E_b[] = {1., 15.055, 48.061, 99.304, 0.012, 8.650};

double arr_M_a[] = {1., -1.465, 1.260, 0.262};
double arr_M_b[] = {1., 9.627, 0.000, 0.000, 11.179, 13.245};

double puc_E_a[] = {1., -0.299};
double puc_E_b[] = {1., 11.11, 14.11, 15.7};

double puc_M_a[] = {1., 0.081};
double puc_M_b[] = {1., 11.15, 18.45, 5.31};

double pol(double p[], unsigned int N, double tau)
{
	double s = 0.;
	double f = 1.;
	for (unsigned int i = 0; i < N; i++, f *= tau)
		s += p[i] * f;

	return s;
}

//----------------------------------------------------------------------------------------------------

double CoulombInterference::FF_dipole(double t) const
{
	double f = (di_la_sq / (di_la_sq - t));
	return f*f;
}

//----------------------------------------------------------------------------------------------------

double CoulombInterference::FF_e(double t) const
{
	/// t is negative
	
	double tau = -t/4./cnts->M_sq;

	switch (ffType)
	{
		case ffNone:
			return 1.;

		case ffDipole:
			return FF_dipole(t);

		case ffHofstadter:
			return FF_dipole(t) * (1. - tau*cnts->kappa);

		case ffBorkowski:
			return bor_E_c[0] / (bor_E_msq[0] - t) 
				+ bor_E_c[1] / (bor_E_msq[1] - t) 
				+ bor_E_c[2] / (bor_E_msq[2] - t) 
				+ bor_E_c[3] / (bor_E_msq[3] - t);
		
		case ffKelly:
			return pol(kel_E_a, 2, tau) / pol(kel_E_b, 4, tau);

		case ffArrington:
			return pol(arr_E_a, 4, tau) / pol(arr_E_b, 6, tau);

		case ffPuckett:
			return pol(puc_E_a, 2, tau) / pol(puc_E_b, 4, tau);
		
		case ffPuckettEl:
			return pol(puc_E_a, 2, tau) / pol(puc_E_b, 4, tau);

		default:
			return 0.;
	}
}

//----------------------------------------------------------------------------------------------------

double CoulombInterference::FF_m(double t) const
{
	/// t is negative
	
	double tau = -t/4./cnts->M_sq;

	switch (ffType)
	{
		case ffNone:
			return 1.;
		
		case ffDipole:
			return FF_dipole(t);
		
		case ffHofstadter:
			return FF_dipole(t) * (1. + cnts->kappa);

		case ffBorkowski:
			return (1. + cnts->kappa) *
				(bor_M_c[0] / (bor_M_msq[0] - t) 
				+ bor_M_c[1] / (bor_M_msq[1] - t) 
				+ bor_M_c[2] / (bor_M_msq[2] - t) 
				+ bor_M_c[3] / (bor_M_msq[3] - t));
		
		case ffKelly:
			return (1. + cnts->kappa) * pol(kel_M_a, 2, tau) / pol(kel_M_b, 4, tau);

		case ffArrington:
			return (1. + cnts->kappa) * pol(arr_M_a, 4, tau) / pol(arr_M_b, 6, tau);

		case ffPuckett:
			return (1. + cnts->kappa) * pol(puc_M_a, 2, tau) / pol(puc_M_b, 4, tau);
		
		case ffPuckettEl:
			return 0.;

		default:
			return 0.;
	}
}

//----------------------------------------------------------------------------------------------------

double CoulombInterference::FF_sq(double t) const
{
	double eff = FF_e(t);

	if (ffType == ffPuckettEl)
		return eff*eff;

	double mff = FF_m(t);
	double tau = -t/4./cnts->M_sq;

	return (eff*eff + tau*mff*mff) / (1. + tau);	
}

//----------------------------------------------------------------------------------------------------

string CoulombInterference::GetFFName() const
{
	switch (ffType)
	{
		case ffNone: return "none";
		case ffDipole: return "dipole";
		case ffHofstadter: return "Hofstadter";
		case ffBorkowski: return "Borkowski";
		case ffKelly: return "Kelly";
		case ffArrington: return "Arrington";
		case ffPuckett: return "Puckett";
		case ffPuckettEl: return "Puckett electric only";
		default: return "unknown";
	}
}

//------------------------------------- COULOMB AMPLITUDE ------------------------------------------

TComplex CoulombInterference::Amp_pure(double t) const
{
	switch (cnts->pMode)
	{
		case Constants::mPP:
			return + cnts->alpha * cnts->s * FF_sq(t) / t;
		case Constants::mAPP:
			return - cnts->alpha * cnts->s * FF_sq(t) / t;
		default:
			return 0.;
	}
}

//--------------------------------------------------------------------------------------------------

double CoulombInterference::A_integrand(double tt, double *par, const void *vobj)
{
	CoulombInterference *obj = (CoulombInterference *) vobj;
	const double &t = par[0];

	return log(tt / t) * obj->FF_sq_prime(tt);
}

//--------------------------------------------------------------------------------------------------

double CoulombInterference::A_term(double t) const
{
	double par[] = { t };
	return RealIntegrate(A_integrand, par, this, t-T, 0., 0., precision, integ_workspace_size,
		integ_workspace, "CoulombInterference::A_term");
}

//--------------------------------------------------------------------------------------------------

double CoulombInterference::I_integrand(double phi, double *par, const void *vobj)
{
	CoulombInterference *obj = (CoulombInterference *) vobj;

	const double &tp = par[0];
	const double &t = par[1];
	
	// t2 ... t''
	double t2 = tp + t + 2. * sqrt(tp * t) * cos(phi);
	return obj->FF_sq(t2) / t2;
}

//--------------------------------------------------------------------------------------------------

double CoulombInterference::I_function_integration(double t, double tp) const
{
	double par[] = { tp, t };

	return RealIntegrate(I_integrand, par, this, 0., 2.*cnts->pi, 0., precision, integ_workspace_size,
		integ_workspace2, "CoulombInterference::I_integration");
}

//--------------------------------------------------------------------------------------------------

void CoulombInterference::InitIFunctionInterpolator(double mt_max, unsigned int n_grid_values)
{
	iFunctionInterpolator = new IFunctionInterpolator(this, mt_max, n_grid_values);
	useIFunctionInterpolator = true;
}

//--------------------------------------------------------------------------------------------------

double CoulombInterference::I_function(double t, double tp) const
{
	if (useIFunctionInterpolator)
		return iFunctionInterpolator->Eval(t, tp);
	else
		return I_function_integration(t, tp);
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::B_integrand(double tp, double *par, const void *vobj)
{
	CoulombInterference *obj = (CoulombInterference *) vobj;
	
	const double &t = par[0];
	TComplex T_hadron_t(par[1], par[2]);

	double I = obj->I_function(t, tp);
	
	TComplex a = model->Amp(tp) / T_hadron_t;
	return (a - 1.) * I	/ 2. / cnts->pi;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::B_term(double t) const
{
	/**
	Function B_integrand(t', t) has problems at point t' = t. It is not defined there (left and right limits are different),
	it is not continuous at the point. To avoid problems, we cut out a small interval (t-tau, t+tau), tau > 0 from the
	the integration region (t_min, 0). The second note concerns exponential fall off of B_integrand as one goes with t' away
	from t. Thus one can take (t-T, 0) instead of (t_min, 0). Of course, T must be suffciently large.
	
	In fact, one must be careful with lower bound t+tau, since it must be less than 0. Otherwise contribution from the 
	region (t+tau, 0) isn't present.
	**/

	TComplex amp_t = model->Amp(t);
	double par[] = { t, amp_t.Re(), amp_t.Im() }; 

	TComplex I = ComplexIntegrate(B_integrand, par, this, t - T, t - tau, 1E-7, precision,
		integ_workspace_size, integ_workspace, "CoulombInterference::B_term/ht");
	if (t + tau < 0.)
		I += ComplexIntegrate(B_integrand, par, this, t + tau, 0., 1E-7, precision,
			integ_workspace_size, integ_workspace, "CoulombInterference::B_term/lt");

	return I;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::B_cahn_integrand_phi(double phi, double *par, const void* /*vobj*/)
{
	const double &tp = par[0];
	const double &t = par[1];
	TComplex T_hadron_t(par[2], par[3]);
	
	double tpp = tp + t + 2. * sqrt(tp * t) * cos(phi);

	TComplex ret = model->Amp(tpp) / T_hadron_t - 1.;

	// The above calculation of ret suffers from round-off errors for constant phase:
	// ret value fluctuates between zero and some very small number.
	// The code below is a work around.
	if (fabs(ret.Im()) < 1E-11)
		ret = TComplex(ret.Re(), 0.);

	return ret;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::B_cahn_integrand_t(double tp, double *par, const void *vobj)
{
	CoulombInterference *obj = (CoulombInterference *) vobj;
	
	const double &t = par[0];

	double ppar[] = { tp, t, par[1], par[2] };
	TComplex I = ComplexIntegrate(B_cahn_integrand_phi, ppar, vobj, 0., 2.*cnts->pi, 0., obj->precision,
		obj->integ_workspace_size, obj->integ_workspace2, "B_cahn_integrand_t");
	
	return obj->FF_sq(tp) / tp / 2. / cnts->pi * I;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::B_term_cahn(double t) const
{
	TComplex amp_t = model->Amp(t);
	double par[] = { t, amp_t.Re(), amp_t.Im() }; 

	TComplex I = ComplexIntegrate(B_cahn_integrand_t, par, this, t - T, t - tau, 1E-7, precision,
		integ_workspace_size, integ_workspace, "CoulombInterference::B_term_cahn/ht");
	if (t + tau < 0.)
		I += ComplexIntegrate(B_cahn_integrand_t, par, this, t + tau, 0., 1E-7, precision,
			integ_workspace_size, integ_workspace, "CoulombInterference::B_term_cahn/lt");

	return I;
}

//--------------------------------------------------------------------------------------------------

double CoulombInterference::C_term(double t) const
{
	return FF_sq(t-T) * log(t/(t-T));
	//return FF_sq(cnts->t_min) * log(t/cnts->t_min);
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Phi_WY(double t) const
{
	// TODO: implement Eq. (23) in [1], but with "-" in front of eta in the second term
	TComplex corr = (cnts->pMode == cnts->mPP) ? -cnts->alpha : +cnts->alpha;
	corr *= + B_term(t) + C_term(t);
	return corr;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Phi_SWY(double t) const
{
	// step for discrete derivative
	double zero = -1E-7;
	double ep = -1E-5;

	// diffractive slope (at t=0)
	double B = log(model->Amp(zero+ep).Rho2() / model->Amp(zero).Rho2()) / ep;

	double phi = cnts->gamma + log(- B * t / 2.);
	if (cnts->pMode == cnts->mPP)
		phi = -phi;

	return cnts->alpha * phi;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Psi_KL(double t) const
{
	TComplex corr = (cnts->pMode == cnts->mPP) ? -cnts->alpha : +cnts->alpha;
	corr *= A_term(t) - B_term(t) - C_term(t);
	return corr;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Psi_Cahn(double t) const
{
	TComplex corr = (cnts->pMode == cnts->mPP) ? -cnts->alpha : +cnts->alpha;
	corr *= A_term(t) - B_term_cahn(t);
	return corr;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Phase(double t) const
{
	switch (mode)
	{
		case mPC: return 0.;
		case mPH: return 0.;
		case mWY: return -Phi_WY(t);
		case mSWY: return -Phi_SWY(t);
		case mKL: return Psi_KL(t);
		case mCahn: return Psi_Cahn(t);
		default: return 0.;
	};
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Amp_WY(double t) const
{
	return Amp_pure(t) * TComplex::Exp(i*Phi_WY(t)) + model->Amp(t);
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Amp_SWY(double t) const
{
	return Amp_pure(t) * TComplex::Exp(i*Phi_SWY(t)) + model->Amp(t);
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Amp_KL(double t) const
{
	return Amp_pure(t) + model->Amp(t) * TComplex::Exp(i*Psi_KL(t));
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Amp_Cahn(double t) const
{
	return Amp_pure(t) + model->Amp(t) * TComplex::Exp(i*Psi_Cahn(t));
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Amp(double t) const
{
	switch (mode)
	{
		case mPC:
			return Amp_pure(t);
		case mPH:
			return model->Amp(t);
		case mWY:
			return Amp_WY(t);
		case mSWY:
			return Amp_SWY(t);
		case mKL:
			return Amp_KL(t);
		case mCahn:
			return Amp_Cahn(t);
		default:
			return 0.;
	};
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::R(double t) const
{
	double KL = Amp_KL(t).Rho2();
	double SWY = Amp_SWY(t).Rho2();
	return (KL - SWY) / KL;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::R_with_cutoff(double t, double cutoff) const
{
	double KL = Amp_KL(t).Rho2();
	double ref = (t > cutoff) ? Amp_WY(t).Rho2() : model->Amp(t).Rho2() ;
	return (KL - ref) / KL;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::Z(double t) const
{
	double KL = Amp_KL(t).Rho2();
	double PH = model->Amp(t).Rho2();
	double PC = Amp_pure(t).Rho2();
	return (KL - PH - PC) / KL;
}

//--------------------------------------------------------------------------------------------------

TComplex CoulombInterference::C(double t) const
{
	double KL = Amp_KL(t).Rho2();
	double PH = model->Amp(t).Rho2();
	return (KL - PH) / PH;
}

} // namespace
