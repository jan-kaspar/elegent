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

#include "interface/PPPModel.h"
#include "interface/Constants.h"

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

PPPModel::PPPModel()
{
	fullLabel.name = "Petrov et al."; shortLabel.name = "petrov";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

PPPModel::~PPPModel()
{
	if (integ_workspace_initialized)
	{
		gsl_integration_workspace_free(integ_workspace);
	}
}

//----------------------------------------------------------------------------------------------------

void PPPModel::Configure(PPPModel::VariantType _v)
{
	variant = _v;

	if (variant == v2P)
	{
		fullLabel.variant = "2P"; shortLabel.variant = "2p";
	}
	
	if (variant == v3P)
	{
		fullLabel.variant = "3P"; shortLabel.variant = "3p";
	}

	if (variant != v2P && variant != v3P)
		printf("ERROR in PPPModel::Init > Unknown variant %u.\n", variant);

	// set labels
	fullLabel.version = "Eur. Phys. J. C23 (2002) 135-143"; shortLabel.version = "02";
	fullLabel.mode = ""; shortLabel.mode = "";
}

//----------------------------------------------------------------------------------------------------

void PPPModel::SetTrajectory(Trajectory &t, double D, double c, double ap, double r2, double s0)
{
	t.D = D; t.c = c; t.ap = ap; t.r2 = r2;
	
	// now pre-compute rho2 and gamma factor
	// physics must be already intialized !!!
	t.rho2 = 4. * t.ap * log(cnts->s/s0) + t.r2;
	t.gamma = t.c / s0 * TComplex::Power(-i * cnts->s/s0, t.D);
}

//----------------------------------------------------------------------------------------------------

void PPPModel::Init()
{
	// physics parameters
	s0 = 1.;	// in GeV^2
	
	//						Delta		c			a'		r^2
	//						1			1			GeV^2	GeV^2
	if (variant == v2P)
	{
		// parameters from Table 1 in [1]
		SetTrajectory(pom1, 0.08590,	53.18,		0.360,	9.595,	s0);
		SetTrajectory(pom2, 0.14437,	6.87,		0.082,	4.765,	s0);
		SetTrajectory(oder, -0.2707,	1.8134,		0.029,	1.159,	s0);
		SetTrajectory(regf, -0.3100,	188.51,		0.84,	41.424, s0);
		SetTrajectory(rego, -0.5300,	-171.36,	0.93,	2.621,	s0);
	}
	
	if (variant == v3P)
	{
		// parameters from Table 2 in [1]
		SetTrajectory(pom1, 0.0578,		53.007,	 	0.5596, 6.3096, s0);
		SetTrajectory(pom2, 0.1669,		9.6762,		0.2733,	3.1097, s0);
		SetTrajectory(pom3, 0.2032,		1.6654,		0.0937,	2.4771, s0);
		SetTrajectory(oder, 0.1920,		0.0166,		0.048,	0.1398, s0);
		SetTrajectory(regf, -0.31,		191.69,		0.84,	31.593, s0);
		SetTrajectory(rego, -0.53,		-174.18,	0.93,	7.467,	s0);
	}
	
	// integration parameters
	precision = 1E-2;
	upper_bound = 40.;

	// prepare integration workspace
	if (!integ_workspace_initialized)
	{
		integ_workspace_size = 100;
		integ_workspace = gsl_integration_workspace_alloc(integ_workspace_size);
		integ_workspace_initialized = true;
	}
}

//----------------------------------------------------------------------------------------------------

void PPPModel::Print() const
{
	printf(">> PPPModel::Print\n");
	printf("\t%s\n", CompileFullLabel().c_str());
	printf("\tvariant: %u\n", variant);

	printf("\tpom1: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", pom1.D, pom1.c, pom1.ap, pom1.r2);
	printf("\tpom2: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", pom2.D, pom2.c, pom2.ap, pom2.r2);
	if (variant == v3P)
		printf("\tpom3: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", pom3.D, pom3.c, pom3.ap, pom3.r2);
	printf("\toder: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", oder.D, oder.c, oder.ap, oder.r2);
	printf("\tregf: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", regf.D, regf.c, regf.ap, regf.r2);
	printf("\trego: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", rego.D, rego.c, rego.ap, rego.r2);
	printf("\tupper_bound=%.1f, precision=%.1E\n", upper_bound, precision);
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::Delta(const Trajectory &traj, double b)
{
	/// delta+-(s, b) according Eq. (11) without the leading i and Eq. (12)

	// NOTE: the s-dependence is calculated during initialization, therefore any cnts->s change
	// after Init() will have no effect

	return traj.gamma * exp(- b * b / traj.rho2) / (4. * cnts->pi * traj.rho2);
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::prf0(double b) const
{
	/// delta as from Eq. (9) from [1]
	TComplex delta = i*Delta(pom1, b) + i*Delta(pom2, b) + i*Delta(regf, b);

	if (cnts->pMode == cnts->mPP)
		delta += Delta(oder, b) + Delta(rego, b);
	else
		delta -= Delta(oder, b) + Delta(rego, b);

	if (variant == v3P)
		delta += i*Delta(pom3, b);

	/// scattering amplitude according to Eq. (1) in [1]
	return (TComplex::Exp(2.*i*delta) - 1.) / 2. / i;
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::Prf(double b) const
{
	return prf0(b / cnts->hbarc); 
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::prf_J0(double b, double *par, const void *vobj)
{
	const PPPModel *obj = (PPPModel *) vobj;
	const double &t = par[0];

	return obj->prf0(b) * b * TMath::BesselJ0(b * sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::Amp(double t) const
{
	double par[] = { t };
	TComplex I = ComplexIntegrate(prf_J0, par, this, 0., upper_bound, precision, integ_workspace_size,
		integ_workspace, "PPPModel::Amp");

	return 2.*cnts->p_cms*cnts->sqrt_s * I;
}
