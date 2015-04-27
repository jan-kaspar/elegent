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

#include "interface/DLModel.h"

#include "interface/Constants.h"
#include "interface/Math.h"

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

DLModel::DLModel()
{
	fullLabel.name = "Donnachie et al."; shortLabel.name = "dl";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

DLModel::~DLModel()
{
	if (integ_workspace_initialized)
	{
		gsl_integration_workspace_free(integ_workspace_t);
	}
}

//----------------------------------------------------------------------------------------------------

void DLModel::Configure()
{
	fullLabel.variant = "";
	fullLabel.version = "Phys. Lett. B727 (2013) 500-505";
	fullLabel.mode = "";
	
	shortLabel.variant = "";
	shortLabel.version = "13";
	shortLabel.mode = "";
}

//----------------------------------------------------------------------------------------------------

void DLModel::Init()
{
	// Eq. (5) in [1]
	ep_P = 0.110;
	ep_pl = -0.327;
	ep_mi = -0.505;

	X_P = 339.;
	X_pl = 212.;
	X_mi = 104.;

	al_P_p = 0.165;		// GeV^-2

	A = 0.682;
	a = 7.854;			// GeV^-2
	b = 2.470;			// GeV^-2

	C = 0.0406;
	t0 = 4.230;			// GeV^2

	// above Eq. (1a) in [1]
	al_pl_p = 0.8;		// GeV^-2
	al_mi_p = 0.92;		// GeV^-2

	// nowhere in [1]
	//nu = 1.;
	nu = cnts->s / 2.;	// TODO: relation guessed

	// integration parameters
	upper_bound_t = 80.; precision_t = 1E-3;

	// prepare integration workspace
	if (!integ_workspace_initialized)
	{
		integ_workspace_size_t = 100;
		integ_workspace_t = gsl_integration_workspace_alloc(integ_workspace_size_t);

		integ_workspace_initialized = true;
	}
}

//----------------------------------------------------------------------------------------------------

void DLModel::Print() const
{
	printf(">> DLModel::Print\n");
	printf("\t%s\n", CompileFullLabel().c_str());

	// TODO

	printf("\n");
	printf("\tintegration parameters:\n");
	printf("\t\tt: upper bound = %.1E, precision = %.1E\n", upper_bound_t, precision_t);
}

//----------------------------------------------------------------------------------------------------

TComplex DLModel::Prf(double b) const
{
	// TODO
	return 0. * b;
}

//----------------------------------------------------------------------------------------------------

double DLModel::F(double t) const
{
	return A * exp(a*t) + (1. - A) * exp(b*t);
}

//----------------------------------------------------------------------------------------------------

TComplex DLModel::A_single(double t) const
{
	double F_P = F(t), F_pl = F_P, F_mi = F_P;

	// Eq. (1a) in [1]
	double al_P = 1. + ep_P + al_P_p*t;
	double al_pl = 1. + ep_pl + al_pl_p*t;
	double al_mi = 1. + ep_mi + al_mi_p*t;

	TComplex A_P = X_P * F_P / 2. / nu * TComplex::Exp(-i*cnts->pi/2.*al_P) * pow(2.*nu*al_P_p, al_P); 
	TComplex A_pl = X_pl * F_pl / 2. / nu * TComplex::Exp(-i*cnts->pi/2.*al_pl) * pow(2.*nu*al_pl_p, al_pl); 
	TComplex A_mi = X_mi * F_mi / 2. / nu * TComplex::Exp(-i*cnts->pi/2.*al_mi) * pow(2.*nu*al_mi_p, al_mi); 

	double sgn = (cnts->pMode = cnts->mPP) ? -1. : +1.;

	return -A_P - A_pl + sgn * i * A_mi;
}

//----------------------------------------------------------------------------------------------------

TComplex DLModel::A_PP(double t) const
{
	double al_PP = 1. + 2.*ep_P + al_P_p/2. * t;

	TComplex L = log(2.*nu*al_P_p) - i*cnts->pi/2.;

	TComplex bracket = A*A/(a + al_P_p*L) * exp(a*t/2.) + (1.-A)*(1.-A)/(b + al_P_p*L) * exp(b*t/2.);

	TComplex orig = X_P*X_P/32./cnts->pi * TComplex::Exp(-i*cnts->pi/2.*al_PP) * pow(2.*nu*al_P_p, al_PP) * bracket;

	return orig / (2.*nu) / 6.;	// TODO: factor needed to reproduce curves
}

//----------------------------------------------------------------------------------------------------

TComplex DLModel::A_ggg(double t) const
{
	if (-t > t0)
		return C * t0*t0*t0 / (t*t*t*t);
	else
		return C / t0 * exp(2. * (1. - t*t/t0/t0));
}

//----------------------------------------------------------------------------------------------------

TComplex DLModel::Amp(double t) const
{
	TComplex A = A_single(t) + A_PP(t) + A_ggg(t);

	// normalisation given below Eq. (1b) in [1]
	return A * cnts->p_cms * cnts->sqrt_s / 4. / cnts->pi;
}
