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

#include "interface/FerreiraModel.h"

#include "interface/Constants.h"
#include "interface/Math.h"

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

FerreiraModel::FerreiraModel()
{
	fullLabel.name = "Ferreira et al."; shortLabel.name = "ferreira";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

FerreiraModel::~FerreiraModel()
{
	if (integ_workspace_initialized)
		gsl_integration_workspace_free(integ_workspace_t);
}

//----------------------------------------------------------------------------------------------------

void FerreiraModel::Configure()
{
	fullLabel.variant = "";
	fullLabel.version = "arXiv:1408.1599v1";
	fullLabel.mode = "";
	
	shortLabel.variant = "";
	shortLabel.version = "14";
	shortLabel.mode = "";
}

//----------------------------------------------------------------------------------------------------

void FerreiraModel::Init()
{
	// page 2 in [1]
	a0 = 1.39;	// GeV^-2

	double W = cnts->sqrt_s / 1E3;	// sqrt(s) in TeV
	double L = log(W / 30.4469);

	// Eq. (28) to (35) in [1]
	al_i = 11.0935 + 1.35479 * log(W);
	al_r = 0.208528 + 0.0419028 * log(W);

	be_i = 4.44606586 + 0.3208411 * L + 0.0613381 * sqrt(0.5 + L*L);
	be_r = 1.1506 + 0.12584 * log(W) + 0.017002 * log(W) * log(W);

	ga_i = 10.025 + 0.79097 * log(W) + 0.088 * log(W) * log(W);
	ga_r = 10.401 + 1.4408 * log(W) + 0.16659 * log(W) * log(W);

	la_i = 14.02008 + 3.23842 * log(W) + 0.444594 * log(W) * log(W);
	la_r = 3.31949 + 0.743706 * log(W);

	N = sqrt(cnts->s * cnts->p_cms * cnts->p_cms / cnts->pi);

	// integration parameters
	precision_t = 1E-4;
	upper_bound_t = -50.;

	if (!integ_workspace_initialized)
	{
		integ_workspace_size_t = 100;
		integ_workspace_t = gsl_integration_workspace_alloc(integ_workspace_size_t);
		integ_workspace_initialized = true;
	}
}

//----------------------------------------------------------------------------------------------------

void FerreiraModel::Print() const
{
	printf(">> FerreiraModel::Print\n");
	printf("\t%s\n", CompileFullLabel().c_str());
	printf("\ta0 = %E GeV^-2\n", a0);
	printf("\tal_i = %E, al_r = %E\n", al_i, al_r);
	printf("\tbe_i = %E, be_r = %E\n", be_i, be_r);
	printf("\tga_i = %E, ga_r = %E\n", ga_i, ga_r);
	printf("\tla_i = %E, la_r = %E\n", la_i, la_r);
}

//----------------------------------------------------------------------------------------------------

double FerreiraModel::Psi(double ga, double t) const
{
	double mt = -t;
	double s1 = sqrt(1. + a0 * mt);
	double s4 = sqrt(4. + a0 * mt);

	return 2. * exp(ga) * (exp(-ga*s1)/s1 - exp(ga) * exp(-ga*s4)/s4);
}

//----------------------------------------------------------------------------------------------------

double FerreiraModel::R_ggg(double t) const
{
	if (fabs(t) < 1E-10)
		return 0.;

	double t2 = t * t;
	double t4 = t2 * t2;
	double r = 0.45/t4 * (1. - exp(-0.005*t4)) * (1. - exp(-0.1*t2));

	return (cnts->pMode == Constants::mPP) ? +r : -r;
}

//----------------------------------------------------------------------------------------------------

TComplex FerreiraModel::Amp(double t) const
{
	double mt = -t;

	// Eq. (18) in [1]
	double re = al_r * exp(-be_r * mt) + la_r * Psi(ga_r, t);
	double im = al_i * exp(-be_i * mt) + la_i * Psi(ga_i, t);

	// Eq. (21) in [1]
	re += R_ggg(t);

	return N * TComplex(re, im);
}

//----------------------------------------------------------------------------------------------------

TComplex FerreiraModel::Amp_J0(double t, double *par, const void *vobj)
{
	const FerreiraModel *obj = (FerreiraModel *) vobj;
	const double &b = par[0];	// impact parameter in GeV^-1

	return obj->Amp(t) * TMath::BesselJ0(b * sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex FerreiraModel::Prf(double b_fm) const
{
	double b = b_fm / cnts->hbarc;	// b in GeV^-1
	double par[] = { b };

	TComplex I = ComplexIntegrate(Amp_J0, par, this, upper_bound_t, 0., 0., precision_t,
		integ_workspace_size_t, integ_workspace_t, "FerreiraModel::Prf");

	return I / 4. / cnts->p_cms / cnts->sqrt_s;
}
