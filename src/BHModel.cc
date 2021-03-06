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

#include "interface/BHModel.h"
#include "interface/Constants.h"

using namespace Elegent;

//----------------------------------------------------------------------------------------------------

BHModel::BHModel()
{
	fullLabel.name = "Block et al."; shortLabel.name = "block";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

BHModel::~BHModel()
{
	if (integ_workspace_initialized)
	{
		gsl_integration_workspace_free(integ_workspace);
	}
}

//----------------------------------------------------------------------------------------------------

void BHModel::Configure()
{
	// set labels
	fullLabel.variant = ""; shortLabel.variant = "";
	fullLabel.version = "Phys. Rept. 436 (2006) 71-215"; shortLabel.version = "06";
	fullLabel.mode = ""; shortLabel.mode = "";
}

//----------------------------------------------------------------------------------------------------

void BHModel::Init()
{
	// parameters from [2]: Table 5, N_g from above Eq. (463), a's and b's from A.1.2

	// common parameters
	m0 = 0.6;
	s0 = 9.53;
	al_s = 0.5;
	Sigma_gg = 9. * cnts->pi * al_s*al_s / m0/m0;

	// sigma_gg
	Cp_gg = 0.00103;
	epsilon = 0.05;

	// Inconsitency in the definition of Ng:
	// [1] : Ng = (6.-epsilon)*(5.-epsilon)*(4.-epsilon)*(3.-epsilon)*(2.-epsilon)*(1.-epsilon) /5./4./3./2. / 2.;	// = 2.649
	// [2] : Ng = 3./2. * (5.-epsilon)*(4.-epsilon)*(3.-epsilon)*(2.-epsilon)*(1.-epsilon) /5./4./3./2.;	// = 1.335
	// Confirmed by F. Halzen: the Ng must be such that "average x carried by the gluons is 1/2",
	//	which corresponds the definition in [1]
	Ng = (6.-epsilon)*(5.-epsilon)*(4.-epsilon)*(3.-epsilon)*(2.-epsilon)*(1.-epsilon) /5./4./3./2. / 2.;
		
	a.push_back(-41.1);
	a.push_back(-487.5);
	a.push_back(-600.);
	a.push_back(600.);
	a.push_back(487.5);
	a.push_back(41.1);

	b.push_back(-9.);
	b.push_back(-225.);
	b.push_back(-900.);
	b.push_back(-900.);
	b.push_back(-225.);
	b.push_back(-9.);

	// sigma_qg
	C_qg_log = 0.166;

	// sigma_qq
	C = 5.36;
	C_even_regge = 29.7;

	// sigma_odd
	C_odd = 10.3;
	
	// mu's
	mu_gg = 0.73;
	mu_odd = 0.53;
	mu_qq = 0.89;

	// precompute mu_qg
	mu_qg = sqrt(mu_qq*mu_gg);
		
	// precompute sigma_gg, Eq. (B5) in [1]: below the Cp_gg stands for C'_gg
	sigma_gg = Cp_gg * Sigma_gg * Ng*Ng * Sum(cnts->s);

	// precompute sigma_qq, Eq. (B9) in [1]
	sigma_qq = Sigma_gg * (C + C_even_regge * m0/sqrt(cnts->s) * TComplex::Exp(i * cnts->pi / 4.));

	// precompute sigma_qg, Eq. (B10) in [1]
	sigma_qg = Sigma_gg * C_qg_log * TComplex(log(cnts->s/s0), -cnts->pi/2.);
	
	// precompute sigma_odd, Eq. (B12) in [1] - the factor in front of W(b, mu_odd)
	// plus additional factor (-i) to match the normalization used in chi_without_i
	sigma_odd = -i * C_odd * Sigma_gg * m0 / cnts->sqrt_s * TComplex::Exp(i * cnts->pi / 4.);

	// integration settings
	upper_bound = 50.;
	precision = 1E-3;

	// prepare integration workspace
	if (!integ_workspace_initialized)
	{
		integ_workspace_size = 100;
		integ_workspace = gsl_integration_workspace_alloc(integ_workspace_size);
		integ_workspace_initialized = true;
	}
}

//----------------------------------------------------------------------------------------------------

void BHModel::Print() const
{
	printf(">> BHModel::Print\n");
	printf("\t%s\n", CompileFullLabel().c_str());

	printf("\tcommon:\n");
	printf("\t\ts0 = %E\n", s0);
	printf("\t\tm0 = %E\n", m0);
	printf("\t\tSigma_gg = %E\n", Sigma_gg);

	printf("\tsigma_gg:\n");
	printf("\t\tCp_gg = %E\n", Cp_gg);
	printf("\t\tNg = %E\n", Ng);
	printf("\t\tepsilon = %E\n", epsilon);
	printf("\t\ta0 = %E, a1 = %E, a2 = %E, a3 = %E, a4 = %E, a5 = %E\n", a[0], a[1], a[2], a[3], a[4], a[5]);
	printf("\t\tb0 = %E, b1 = %E, b2 = %E, b3 = %E, b4 = %E, b5 = %E\n", b[0], b[1], b[2], b[3], b[4], b[5]);
	printf("\t\tsigma_gg: Re = %E, Im = %E\n", sigma_gg.Re(), sigma_gg.Im());

	printf("\tsigma_qg:\n");
	printf("\t\tC_qg_log = %E\n", C_qg_log);
	printf("\t\tsigma_qg: Re = %E, Im = %E\n", sigma_qg.Re(), sigma_qg.Im());

	printf("\tsigma_qq:\n");
	printf("\t\tC = %E\n", C);
	printf("\t\tC_even_regge = %E\n", C_even_regge);
	printf("\t\tsigma_qq: Re = %E, Im = %E\n", sigma_qq.Re(), sigma_qq.Im());

	printf("\tsigma_odd:\n");
	printf("\t\tC_odd = %E\n", C_odd);
	printf("\t\tsigma_odd: Re = %E, Im = %E\n", sigma_odd.Re(), sigma_odd.Im());
	
	printf("\tmu's:\n");
	printf("\t\tmu_gg = %E\n", mu_gg);
	printf("\t\tmu_qq = %E\n", mu_qq);
	printf("\t\tmu_qg = %E\n", mu_qg);
	printf("\t\tmu_odd = %E\n", mu_odd);
	
	printf("\tintegration parameters:\n");
	printf("\t\tupper_bound = %.lf\n", upper_bound);
	printf("\t\tprecision = %.lE\n", precision);
}

//----------------------------------------------------------------------------------------------------

TComplex BHModel::Sum(double s) const
{
	TComplex log_tau0 = log(m0*m0 / s) + i * cnts->pi / 2.;

	TComplex S = 0.;

	for (unsigned int ii = 0; ii <= 5; ii++)
	{
		double i = double(ii);
		double ime = i - epsilon;
		double f = b[ii] / ime;
		double t1 = (a[i] - f) / ime;
		TComplex v = t1 - TComplex::Exp(ime * log_tau0) * (t1 + f * log_tau0);
		S += v;
	}

	return S;
}

//----------------------------------------------------------------------------------------------------
double BHModel::W(double b, double mu) const
{
	// Eq. (B2) in [1]
	double mub = mu * b;

	// evaluates v = (mu b)^3 K_3(mu b); directly or in continuous limit
	double v = 0.;
	if (mub < 1E-10)
		v = 8.;
	else
		v = mub*mub*mub * TMath::BesselK(3, mub);
	return mu*mu * v / 96. / cnts->pi;
}

//----------------------------------------------------------------------------------------------------

TComplex BHModel::chi_without_i(double b) const
{
	// Eqs. (B1) without the leading i factor and Eq. (B12) in [1]
	double odd_sign = (cnts->pMode == cnts->mAPP) ? +1. : -1.;
	return (
			sigma_gg * W(b, mu_gg)
			+ sigma_qg * W(b, mu_qg)
			+ sigma_qq * W(b, mu_qq)
			+ odd_sign * sigma_odd * W(b, mu_odd)
		) / 2.;
}

//----------------------------------------------------------------------------------------------------

TComplex BHModel::prf0(double b) const
{
	return (1. - TComplex::Exp(-chi_without_i(b)));
}

//----------------------------------------------------------------------------------------------------

TComplex BHModel::Prf(double b) const
{
	return prf0(b / cnts->hbarc) * i / 2.;
}

//----------------------------------------------------------------------------------------------------

TComplex BHModel::prf0_J0(double b, double *par, const void *vobj)
{
	const BHModel *obj = (BHModel *) vobj;
	const double &q = par[0];

	return b * TMath::BesselJ0(b*q) * obj->prf0(b);
} 

//----------------------------------------------------------------------------------------------------

TComplex BHModel::Amp(double t) const
{
	// from Eqs. (A11) and (A12)
	double q = sqrt(-t);
	double par[] = { q };
	TComplex I = ComplexIntegrate(prf0_J0, par, this, 0., upper_bound, 0., precision,
		integ_workspace_size, integ_workspace, "BHModel::Amp");
	return i * cnts->p_cms * cnts->sqrt_s * I;
}
