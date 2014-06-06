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

#include "interface/GodizovModel.h"
#include "interface/Constants.h"

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

GodizovModel::GodizovModel()
{
	fullLabel.name = "Godizov"; shortLabel.name = "godizov";
}

//----------------------------------------------------------------------------------------------------

void GodizovModel::Configure()
{
	fullLabel.variant = "";
	fullLabel.version = "arXiv:1404.2851v2";
	fullLabel.mode = "";
	
	shortLabel.variant = "";
	shortLabel.version = "14";
	shortLabel.mode = "";
}

//----------------------------------------------------------------------------------------------------

void GodizovModel::Init()
{
	// TODO: verify - not in the reference paper [1]
	s0 = 1.;	// in GeV^2
	
	// parameters from Table 1 from [1]
	De = 0.111;		// al_P(0) - 1
	ta_a = 0.47;	// GeV^2
	Ga_P0 = 7.43;	
	ta_g = 0.98;	// GeV^2

	// TODO
	gsl_w_size = 1000000;	// TODO: tune
	gsl_w = gsl_integration_workspace_alloc(gsl_w_size);
	// TODO: needs to be freed!!
}

//----------------------------------------------------------------------------------------------------

void GodizovModel::Print() const
{
	printf(">> GodizovModel::Print\n");
	printf("\t%s\n", CompileFullLabel().c_str());
	printf("\tal_p(0) - 1 = %.3f\n", De);
	printf("\ttau_a = %.3f\n", ta_a);
	printf("\tGa_P(0) = %.3f\n", Ga_P0);
	printf("\ttau_g = %.3f\n", ta_g);
	printf("\ts0 = %.3f\n", s0);
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::delta_t(double t) const
{
	/// Eq. (4) in [1]
	double al_P = 1. + De / (1. - t/ta_a);
	double Ga_P = Ga_P0 / (1. - t/ta_g) / (1. - t/ta_g);

	/// delta(s, t) according Eq. (2) in [1]

	return (i + tan(cnts->pi * De / 2.)) * Ga_P*Ga_P * pow(cnts->s / s0, al_P);
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::delta_t_J0(double t, void *vpa)
{
	void **vp = (void **) vpa;
	const GodizovModel *obj = (GodizovModel *) vp[0];

	double *par = (double *) vp[1];
	double b = par[0];
	
	return obj->delta_t(t) * TMath::BesselJ0(b * sqrt(-t));
}

double GodizovModel::delta_t_J0_Re(double t, void *vpa) { return GodizovModel::delta_t_J0(t, vpa).Re(); }
double GodizovModel::delta_t_J0_Im(double t, void *vpa) { return GodizovModel::delta_t_J0(t, vpa).Im(); }

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::delta_b(double b) const
{
	// bottom relation from Eq. (1) in [1]

	// TODO: tune
	double precision_GSL = 1E-4;
	double upper_bound = 60.;	// GeV^-2
	
	double result_re, result_im;

	const void* par[] = { this, &b };

	{
		double error;
		gsl_function F;
	  	F.function = delta_t_J0_Re;
	  	F.params = par;
	
		gsl_integration_qag(&F, -upper_bound, 0., 0., precision_GSL, gsl_w_size, GSL_INTEG_GAUSS61, gsl_w, &result_re, &error);
	}

	{
		double error;
		gsl_function F;
	  	F.function = delta_t_J0_Im;
	  	F.params = par;
	
		gsl_integration_qag(&F, -upper_bound, 0., 0., precision_GSL, gsl_w_size, GSL_INTEG_GAUSS61, gsl_w, &result_im, &error);
	}

	TComplex I(result_re, result_im);

	return I / 16. / cnts->pi / cnts->s;
}
	
//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::prf0(double b) const
{
	//printf(">> GodizovModel::prf0(%E)\n", b);

	return (TComplex::Exp(2. * i * delta_b(b)) - 1.) / (2. * i);
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::Amp(double t) const
{
	//printf(">> GodizovModel::Amp(%E)\n", t);

	// integral from upper part of Eq. (1) in [1]
	// TODO: tune
	double precision_GSL = 1E-4;
	double upper_bound = 50.;	// GeV^-1
	
	double result_re, result_im;

	const void* par[] = { this, &t };

	{
		double error;
		gsl_function F;
	  	F.function = prf_J0_Re;
	  	F.params = par;
	
		gsl_integration_qag(&F, 0., upper_bound, 0., precision_GSL, gsl_w_size, GSL_INTEG_GAUSS61, gsl_w, &result_re, &error);
	}

	{
		double error;
		gsl_function F;
	  	F.function = prf_J0_Im;
	  	F.params = par;
	
		gsl_integration_qag(&F, 0., upper_bound, 0., precision_GSL, gsl_w_size, GSL_INTEG_GAUSS61, gsl_w, &result_im, &error);
	}

	TComplex I(result_re, result_im);

	return 2. * cnts->p_cms * cnts->sqrt_s * I;
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::Prf(double b) const
{
	return prf0(b / cnts->hbarc); 
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::prf_J0(double b, void *vpa)
{
	void **vp = (void **) vpa;
	const GodizovModel *obj = (GodizovModel *) vp[0];

	double *par = (double *) vp[1];
	double t = par[0];

	return obj->prf0(b) * b * TMath::BesselJ0(b * sqrt(-t));
}

double GodizovModel::prf_J0_Re(double b, void *vp) { return GodizovModel::prf_J0(b, vp).Re(); }
double GodizovModel::prf_J0_Im(double b, void *vp) { return GodizovModel::prf_J0(b, vp).Im(); }
