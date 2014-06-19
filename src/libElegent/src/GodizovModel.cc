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
#include "interface/Math.h"

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

GodizovModel::GodizovModel()
{
	fullLabel.name = "Godizov"; shortLabel.name = "godizov";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

GodizovModel::~GodizovModel()
{
	if (integ_workspace_initialized)
	{
		gsl_integration_workspace_free(integ_workspace_b);
		gsl_integration_workspace_free(integ_workspace_t);
	}
}

//----------------------------------------------------------------------------------------------------

void GodizovModel::Configure(bool _presampled)
{
	presampled = _presampled;

	fullLabel.variant = "";
	fullLabel.version = "Phys. Lett. B735 (2014) 57-61";
	fullLabel.mode = "";
	
	shortLabel.variant = "";
	shortLabel.version = "14";
	shortLabel.mode = "";
}

//----------------------------------------------------------------------------------------------------

void GodizovModel::Init()
{
	// neither in [1] nor [2], but confirmed by the author
	s0 = 1.;	// in GeV^2
	
	// parameters from Table 1 from [2]
	De = 0.111;		// al_P(0) - 1
	ta_a = 0.47;	// GeV^2
	Ga_P0 = 7.43;	
	ta_g = 0.98;	// GeV^2

	// integration parameters
	upper_bound_t = 80.; precision_t = 1E-3;
	upper_bound_b = 60.; precision_b = 1E-2;

	// prepare integration workspace
	if (!integ_workspace_initialized)
	{
		integ_workspace_size_t = 100;
		integ_workspace_t = gsl_integration_workspace_alloc(integ_workspace_size_t);

		integ_workspace_size_b = 1000;
		integ_workspace_b = gsl_integration_workspace_alloc(integ_workspace_size_b);

		integ_workspace_initialized = true;
	}

	// pre-sample profile function
	if (presampled)
		Prf0SampleBuild(40001);
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

	printf("\n");
	printf("\tpresampled = %u\n", presampled);
	if (presampled)
		printf("\t\tsample size = %u, db = %.1E\n", prf0_sample_N, prf0_sample_db);
	
	printf("\n");
	printf("\tintegration parameters:\n");
	printf("\t\tt: upper bound = %.1E, precision = %.1E\n", upper_bound_t, precision_t);
	printf("\t\tb: upper bound = %.1E, precision = %.1E\n", upper_bound_b, precision_b);
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::delta_t(double t) const
{
	/// Eq. (5) in [2]
	double al_P = 1. + De / (1. - t/ta_a);
	double Ga_P = Ga_P0 / (1. - t/ta_g) / (1. - t/ta_g);

	/// delta(s, t) according Eq. (3) in [2]
	return (i + tan(cnts->pi * De / 2.)) * Ga_P*Ga_P * pow(cnts->s / s0, al_P);
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::delta_t_J0(double t, double *par, const void *vobj)
{
	const GodizovModel *obj = (GodizovModel *) vobj;
	double b = par[0];

	return obj->delta_t(t) * TMath::BesselJ0(b * sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::delta_b(double b) const
{
	// bottom relation from Eq. (1) in [2]
	double par[] = { b };
	TComplex I = ComplexIntegrate(delta_t_J0, par, this, -upper_bound_t, 0., 0., precision_t,
		integ_workspace_size_t, integ_workspace_t, "GodizovModel::delta_b");
	return I / 16. / cnts->pi / cnts->s;
}
	
//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::prf0(double b) const
{
	return (TComplex::Exp(2. * i * delta_b(b)) - 1.) / (2. * i);
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::Prf(double b) const
{
	return prf0(b / cnts->hbarc); 
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::prf_J0(double b, double *par, const void *vobj)
{
	const GodizovModel *obj = (GodizovModel *) vobj;
	double t = par[0];

	TComplex prf0_v = (obj->presampled) ? obj->Prf0SampleEval(b) : obj->prf0(b);
	return prf0_v * b * TMath::BesselJ0(b * sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::Amp(double t) const
{
	// integral from upper part of Eq. (1) in [2]
	double par[] = { t };
	TComplex I = ComplexIntegrate(prf_J0, par, this, 0., upper_bound_b, 0., precision_b,
		integ_workspace_size_b, integ_workspace_b, "GodizovModel::Amp");

	return 2. * cnts->p_cms * cnts->sqrt_s * I;
}
//----------------------------------------------------------------------------------------------------

void GodizovModel::Prf0SampleBuild(unsigned int samples)
{
	printf(">> GodizovModel::Prf0SampleBuild > Building %u samples...\n", samples);

	prf0_sample_re.clear();
	prf0_sample_re.reserve(samples);
	prf0_sample_im.clear();
	prf0_sample_im.reserve(samples);

	double db = upper_bound_b / (samples - 1);
	prf0_sample_db = db;
	prf0_sample_N = samples;

	double b = 0.;
	for (unsigned int i = 0; i < samples; i++, b += db)
	{
		TComplex v = prf0(b);

		prf0_sample_re.push_back(v.Re());
		prf0_sample_im.push_back(v.Im());
	}
}

//----------------------------------------------------------------------------------------------------

TComplex GodizovModel::Prf0SampleEval(double b) const
{
	unsigned int idx = (int)(b / prf0_sample_db);
	
	if (idx + 1 > prf0_sample_N - 1)
		return TComplex(0, 0);

	double f = b/prf0_sample_db - idx;

	return TComplex(
		(prf0_sample_re[idx+1] - prf0_sample_re[idx])*f + prf0_sample_re[idx],
		(prf0_sample_im[idx+1] - prf0_sample_im[idx])*f + prf0_sample_im[idx]
	);
}
