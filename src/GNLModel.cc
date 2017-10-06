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

#include "interface/GNLModel.h"

#include "interface/Constants.h"
#include "interface/Math.h"

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

GNLModel::GNLModel()
{
	fullLabel.name = "Gauron et al."; shortLabel.name = "gnl";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

GNLModel::~GNLModel()
{
	if (integ_workspace_initialized)
	{
		gsl_integration_workspace_free(integ_workspace_t);
	}
}

//----------------------------------------------------------------------------------------------------

void GNLModel::Configure()
{
	fullLabel.variant = "";
	fullLabel.version = "arXiv:hep-ph/0607089v4 (2006)";
	fullLabel.mode = "";
	
	shortLabel.variant = "";
	shortLabel.version = "06";
	shortLabel.mode = "";
}

//----------------------------------------------------------------------------------------------------

void GNLModel::Init()
{
	// Table 1 from [1]
	H_1 = 0.4030;
	b_1_pl = 4.5691;
	H_2 = -3.8616;
	b_2_pl = 7.1798;
	H_3 = 9.2079;
	b_3_pl = 6.0270;
	K_pl = 0.6571;

	O_1 = -0.0690;
	b_1_mi = 8.9526;
	O_2 = 1.4166;
	b_2_mi = 3.4515;
	O_3 = -0.3558;
	b_3_mi = 1.1064;
	K_mi = 0.1267;

	s_P = { 1., 40.43, 4.37, 0.25 };
	s_PP = { 1., -9.2, 1.95, 0 };
	s_O = { 1., -6.07, 5.33, 0.57 };
	s_OP = { 1., 11.83, 1.73, 0 };
	s_R_pl = { 0.48, 38.18, 0.03, 0.88 };
	s_R_mi = { 0.34, 47.09, 33.60, 0.88 };
	s_RP_pl = { -0.56, -1930.1, 0.79, 0. };
	s_RP_mi = { 0.70, 8592.7, 7.33, 0. };

	s_PP.al_p = s_P.al_p / 2.;	// Eq. (16) from [1]
	s_OP.al_p = s_O.al_p * s_P.al_p / (s_O.al_p + s_P.al_p);	// Eq. (31) from [1]
	s_RP_pl.al_p = s_R_pl.al_p * s_P.al_p / (s_R_pl.al_p + s_P.al_p);	// Eq. (22) from [1]
	s_RP_mi.al_p = s_R_mi.al_p * s_P.al_p / (s_R_mi.al_p + s_P.al_p);	// Eq. (37) from [1]

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

void GNLModel::Print() const
{
	printf(">> GNLModel::Print\n");

	// TODO
	
	printf("\n");
	printf("\tintegration parameters:\n");
	printf("\t\tt: upper bound = %.1E, precision = %.1E\n", upper_bound_t, precision_t);
}

//----------------------------------------------------------------------------------------------------

TComplex GNLModel::Prf(double /* b */) const
{
	// TODO
	return 0.; 
}

//----------------------------------------------------------------------------------------------------

TComplex GNLModel::Amp(double t) const
{
	// crossing event part of the amplitude
	const TComplex ln_s_bar = log(cnts->s) - i * M_PI / 2.;
	const TComplex ta_bar = sqrt(-t) * ln_s_bar;

	// TODO: the code below compiles but is wrong - ROOT casts TComplex to double in the arguments of
	// BesselJn
	const TComplex F_H_pl = i * cnts->s * (
			H_1 * ln_s_bar * ln_s_bar * 2. * TMath::BesselJ1(K_pl * ta_bar) / K_pl / ta_bar * exp(b_1_pl * t)
			+ H_2 * ln_s_bar * TMath::BesselJ0(K_pl * ta_bar) * exp(b_2_pl * t)
			+ H_3 * (TMath::BesselJ0(K_pl * ta_bar) - K_pl * ta_bar * TMath::BesselJ1(K_pl * ta_bar)) * exp(b_3_pl * t)
		);

	const double al_P = s_P.al_0 + s_P.al_p * t;
	const TComplex F_P_pl = cnts->s * s_P.C * exp(s_P.be * t) * (i - 1./tan(M_PI/2.*al_P)) * pow(cnts->s, al_P-1.);

	const double al_PP = s_PP.al_0 + s_PP.al_p * t;
	const TComplex F_PP_pl = cnts->s * s_PP.C * exp(s_PP.be * t) * (i*sin(M_PI/2*al_PP) - cos(M_PI/2.*al_PP)) * pow(cnts->s, al_PP-1.) / ln_s_bar;

	const double al_R_pl = s_R_pl.al_0 + s_R_pl.al_p * t;
	const double ga_R_pl = al_R_pl * (al_R_pl+1) * (al_R_pl+2) / s_R_pl.al_0 / (s_R_pl.al_0+1) / (s_R_pl.al_0+2); 
	const TComplex F_R_pl = cnts->s * s_R_pl.C * ga_R_pl * exp(s_R_pl.be) * (i - 1./tan(M_PI/2*al_R_pl)) * pow(cnts->s, al_R_pl-1.);

	const double al_RP_pl = s_RP_pl.al_0 + s_RP_pl.al_p * t;
	const TComplex F_RP_pl = cnts->s *t*t * s_RP_pl.C * exp(s_RP_pl.be * t) * (i*sin(M_PI/2.*al_RP_pl) - cos(M_PI/2.*al_RP_pl))
		* pow(cnts->s, al_RP_pl-1.) / ln_s_bar;

	const TComplex F_pl = F_H_pl + F_P_pl + F_PP_pl + F_R_pl + F_RP_pl;

	// crossing odd part of the amplitude

	const TComplex F_MO_mi = 0.; // TODO

	const TComplex F_O_mi = 0.; // TODO

	const TComplex F_OP_mi = 0.; // TODO

	const TComplex F_R_mi = 0.; // TODO

	const TComplex F_RP_mi = 0.; // TODO

	const TComplex F_mi = F_MO_mi + F_O_mi + F_OP_mi + F_R_mi + F_RP_mi;

	// complete amplitude
	const TComplex F = (cnts->pMode = cnts->mPP) ? (F_pl + F_mi) : (F_pl - F_mi);

	return cnts->p_cms / cnts->sqrt_s / 4. / M_PI * F; 
}
