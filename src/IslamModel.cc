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

#include "interface/IslamModel.h"
#include "interface/Constants.h"

using namespace std;
using namespace Elegent;

//#define DEBUG 1

//----------------------------------------------------------------------------------------------------

IslamModel::IslamModel()
{
	fullLabel.name = "Islam et al."; shortLabel.name = "islam";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

IslamModel::~IslamModel()
{
	if (integ_workspace_initialized)
	{
		gsl_integration_workspace_free(integ_workspace_b);
		gsl_integration_workspace_free(integ_workspace_t);
	}
}

//----------------------------------------------------------------------------------------------------

void IslamModel::Configure(IslamModel::VariantType _v, IslamModel::ModeType _m)
{
	variant = _v;
	mode = _m;
	
	if (variant == vHP)
	{
		fullLabel.variant = "HP"; shortLabel.variant = "hp";
	}
	
	if (variant == vLxG)
	{
		fullLabel.variant = "LxG"; shortLabel.variant = "lxg";
	}

	// set labels
	fullLabel.version = "Int. J. Mod. Phys. A21 (2006) 1-42, Mod. Phys. Lett. A24 (2009) 485-496"; shortLabel.version = "06,09";

	if (mode == mDiff)
		{ fullLabel.mode = "diffraction"; shortLabel.mode = "diff"; }
	if (mode == mCore)
		{ fullLabel.mode = "core"; shortLabel.mode = "core"; }
	if (mode == mQuark)
		{ fullLabel.mode = "quark-quark"; shortLabel.mode = "quark"; }
	if (mode == mDiffCore)
		{ fullLabel.mode = "diffraction+core"; shortLabel.mode = "diff+core"; }
	if (mode == mFull)
		{ fullLabel.mode = "full"; shortLabel.mode = "full"; }
}

//----------------------------------------------------------------------------------------------------

void IslamModel::Init()
{
	// ---------- diffraction amplitude ----------
	// parameters from page 23 of [4] 
	double R0 = 2.77;
   	double R1 = 0.0491;
	double a0 = 0.245;
	double a1 = 0.126;
	R = TComplex(R0 + R1*cnts->ln_s, -R1*cnts->pi/2);
	a = TComplex(a0 + a1*cnts->ln_s, -a1*cnts->pi/2);
	
	double et0 = 0.0844;
	double c0 = 0.0;
	double si = 2.7;
	
	// function g(s) according to Eq. (4.2) in [4]
	TComplex g_s = (1. - CEF(et0, c0, si)) * (1. + TComplex::Exp(-R/a)) / (1. - TComplex::Exp(-R/a));
	Diff_fac = i * cnts->sqrt_s * cnts->p_cms * g_s;
	Diff_fac_profile = g_s * i/2.;
	
	// ---------- absorbtion factor due to diffraction ----------
	
	// parameters from page 23 of [4]
	double la0 = 0.727;
	double d0 = 13.;
	double al = 0.246;
	double hga0 = 1.53;
	double hga1 = 0.;
	double hsi = 1.46;

	if (cnts->pMode == cnts->mPP)
		Abs_fac = cnts->s * CEF(hga0, hga1, hsi) * ( CEF(et0, c0, si) + i*CEF(la0, -d0, al) );

	if (cnts->pMode == cnts->mAPP)
		Abs_fac = cnts->s * CEF(hga0, hga1, hsi) * ( CEF(et0, c0, si) - i*CEF(la0, -d0, al) );
	
	// ---------- core amplitude ----------
	// parameters from page 23 of [4]
	beta = 3.075;
	double m_omega = 0.801;
	m_omega_sq = m_omega * m_omega;
	
	if (cnts->pMode == cnts->mPP)
		Core_fac = -1;

	if (cnts->pMode == cnts->mAPP)
		Core_fac = +1;

	// ---------- quark-quark amplitude ----------
	
	// parameter from page 25 of [4]
	m0sq = 12.;

	// hard pomeron variant
	if (variant == vHP)
	{
		// parameters from page 25 of [4]
		double tgaqq = 0.03;
		omega = 0.15;
		r0 = 2.;

		// the factor (without s) multiplying the 2nd term in the 2nd brackets in Eq. (6.3) in [4]
		Quark_fac = i * tgaqq * TComplex::Power(-i * cnts->s, omega);

		Quark_const = -2. * tgaqq * TComplex::Power(-i * cnts->s, omega);

		// Born term only by default
		qqMaxOrder = 1;
	}
	
	// low-x gluons variant
	if (variant == vLxG)
	{
		double tgagg = 0.0056;	// obtained in private communication with the authors

		// parameters from page 8 of [5]
		lambda = 0.29;
		m_c = 1.67;
		
		// the factor (without is) multiplying the fraction in Eq. (32) in [5]
		cgc_fac = tgagg * TComplex::Power(-i * cnts->s, lambda);

		// Born term only by default
		cgcMaxOrder = 1;
	}
	
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

TComplex IslamModel::CEF(double a, double b, double c)
{
	// crossing-even function:
	//		a + b / (s exp(-i pi/2))^c
	//		= a + b / (-i s)^c = a + b (-i s)^(-c)
	return a + b * TComplex::Power(-i * cnts->s, -c);
}

//----------------------------------------------------------------------------------------------------

void IslamModel::Print() const
{
	printf(">> IslamModel::Print\n");
	printf("\t%s\n", CompileFullLabel().c_str());
	printf("\tdiffraction variables\n");
	double v1 = R.Im(); v1 = -v1 * 2 / cnts->pi;
	double v0 = R.Re(); v0 -= v1 * cnts->ln_s;
	printf("\t\tR0=%f, R1=%f\n", v0, v1);
	v1 = a.Im(); v1 = -v1 * 2 / cnts->pi;
	v0 = a.Re(); v0 -= v1 * cnts->ln_s;
	printf("\t\ta0=%f, a1=%f\n", v0, v1);
	printf("\t\tR: Re=%E, Im=%E\n", R.Re(), R.Im());
	printf("\t\ta: Re=%E, Im=%E\n", a.Re(), a.Im());
	printf("\t\tDiff_fac_profile: Re=%E, Im=%E\n", Diff_fac_profile.Re(), Diff_fac_profile.Im());
	printf("\t\tDiff_fac: Re=%E, Im=%E\n", Diff_fac.Re(), Diff_fac.Im());
	printf("\t\tAbs_fac: Re=%E, Im=%E\n", Abs_fac.Re(), Abs_fac.Im());

	printf("\tcore scattering variables\n");
	printf("\t\tbeta = %E\n", beta);
	printf("\t\tm_omega_sq = %E\n", m_omega_sq);
	printf("\t\tCore_fac = %E\n", Core_fac);

	printf("\tquark-quard scattering variables\n");
	printf("\t\tm0sq = %E\n", m0sq);
	printf("\t\tr0 = %E\n", r0);
	printf("\t\tomega = %E\n", omega);
	printf("\t\tQuark_fac: Re=%E, Im=%E\n", Quark_fac.Re(), Quark_fac.Im());
	printf("\t\tQuark_const: Re=%E, Im=%E\n", Quark_const.Re(), Quark_const.Im());
	printf("\t\tqqMaxOrder = %i\n", qqMaxOrder);
	printf("\t\tlambda = %E\n", lambda);
	printf("\t\tm_c = %E\n", m_c);
	printf("\t\tcgc_fac: Re=%E, Im=%E\n", cgc_fac.Re(), cgc_fac.Im());
	printf("\t\tcgcMaxOrder = %i\n", cgcMaxOrder);

	printf("\tintegration variables\n");
	printf("\t\tprecision_b = %E\n", precision_b);
	printf("\t\tprecision_t = %E\n", precision_t);
	printf("\t\tupper_bound_b = %E\n", upper_bound_b);
	printf("\t\tupper_bound_t = %E\n", upper_bound_t);
}

//-------------------------------------- DIFFRACTION AMPLITUDE ------------------------------------

TComplex IslamModel::GammaD(double b) const
{
	/// b-dependent part of profile function Gamma_D^+ in Eq. (2.7) in [4]
	/// b... impact parameter in	fm
	return 1. / (1. + TComplex::Exp((b - R) / a)) + 1. / (1. + TComplex::Exp((-b - R) / a)) - 1.;
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::GammaD_J0(double b, double *par, const void *vobj)
{
	const IslamModel *obj = (IslamModel *) vobj;
	const double &t = par[0];	

	return obj->GammaD(b) * b * TMath::BesselJ0(b*sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_diff(double t) const
{
#ifdef DEBUG
	printf(">> IslamModel::T_diff\n");
#endif

	/// t < 0
	double par[] = { t };
	return Diff_fac * ComplexIntegrate(GammaD_J0, par, this, 0., upper_bound_b, 0., precision_b,
		integ_workspace_size_b, integ_workspace_b, "IslamModel::T_diff");
}

//----------------------------------------- CORE AMPLITUDE ----------------------------------------

double IslamModel::F_sq(double t)	const
{
	/// formfactor, t < 0
	return beta * sqrt(m_omega_sq - t) * TMath::BesselK1(beta * sqrt(m_omega_sq - t));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_core(double t) const
{
	/// t < 0

#ifdef DEBUG
	printf(">> IslamModel::T_core\n");
#endif

	return Abs_fac * Core_fac * F_sq(t) / (m_omega_sq - t);
}

//---------------------------------------- QUARK AMPLITUDE ----------------------------------------

double IslamModel::I_integral(double qt, double al) const
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

double IslamModel::F_cal_integ(double x, double *par, const void *vobj)
{
	IslamModel *obj = (IslamModel *) vobj;

	const double &qt = par[0];
	const double &n = par[1];
	const double &omega = par[2];
	const double &m0sq = par[3];

	double al_sq = m0sq/4. + cnts->M_sq * x * x;
	double al = sqrt(al_sq);

	return exp((1. + n * omega) * log(x)) / al_sq * obj->I_integral(qt, al);
}

//----------------------------------------------------------------------------------------------------

double IslamModel::F_cal(int n, double qt, double omega, double m0sq) const
{
	double par[] = { qt, double(n), omega, m0sq };
	double I = RealIntegrate(F_cal_integ, par, this, 0., 1., 0., 1E-3, integ_workspace_size_b,
		integ_workspace_b, "IslamModel::F_cal");
	return cnts->M * exp(2.5 * log(m0sq)) / 8. / cnts->pi * I;
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_quark(double t) const
{
	switch (variant)
	{
		case vHP: return T_hp(t);
		case vLxG: return T_lxg(t);
		default:
			printf("ERROR in IslamModel::T_quark > unknown variant %i\n", variant);
			return 0.;
	}
}

//----------------------------------------------------------------------------------------------------

double IslamModel::T_hp_integ(double b, double *par, const void *vobj)
{
	IslamModel *obj = (IslamModel *) vobj;

	const double &q = par[0];
	const double &n = par[1];

	return b * TMath::BesselJ0(b * q) * pow( TMath::BesselK0(b / obj->r0) , n);
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_hp_n(int n, double t) const
{
	double q = sqrt(fabs(t));

	if (n == 1)
		return Quark_fac / (-t + 1./r0/r0);
	
	if (n == 2)
	{
		if (q < 1E-2)
			return i * Quark_fac*Quark_fac * r0*r0*r0/4.;	// limit
		else
			return i * Quark_fac*Quark_fac * 2. * asinh(q * r0 / 2.) / q / sqrt( fabs(t) + 4./r0/r0 );
	}
	
	// general formula
	double par[] = { q, double(n) };
	double I = RealIntegrate(T_hp_integ, par, this, 0., 30., 0., 1E-3, integ_workspace_size_b,
		integ_workspace_b, "IslamModel::T_hp_n");
	return -i / 2. / TMath::Factorial(n) * TComplex::Power(Quark_const, n) * I; // correct -i
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_hp(double t) const
{
#ifdef DEBUG
	printf(">> IslamModel::T_hp\n");
#endif

	/// t < 0
	double qt = sqrt(-t * (-t / cnts->t_min + 1.));

	TComplex sum = 0.;
	for (int j = 1; j <= qqMaxOrder; j++)
	{
		double F = F_cal(j, qt, omega, m0sq);
		sum += F*F * T_hp_n(j, t);
	}

	return Abs_fac * sum;
}

//------------------------------------ CGC AMPLITUDE -------------------------------------------------

double IslamModel::T_lxg_integ(double b, double *par, const void *vobj)
{
	const double &q = par[0];
	const double &n = par[1];

	IslamModel *obj = (IslamModel *) vobj;
	const double &m_c = obj->m_c;

	return b * TMath::BesselJ0(b * q) * pow(exp(-b * m_c) * m_c*m_c * (1. + b*m_c) / 3., n);
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_lxg_n(int n, double t) const
{
	double q = sqrt(fabs(t));

	if (n == 1)
		return i * cgc_fac / pow(1. - t/m_c/m_c, 2.5);
	
	// general formula
	double par[] = { q, double(n) };
	double I = RealIntegrate(T_lxg_integ, par, this, 0., 30., 0., 1E-3, integ_workspace_size_b,
		integ_workspace_b, "IslamModel::T_lxg_n");
	return i * pow(-2., n - 1) / TMath::Factorial(n) * TComplex::Power(cgc_fac, n) * I;
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_lxg(double t) const
{
#ifdef DEBUG
	printf(">> IslamModel::T_lxg\n");
#endif

	/// t < 0
	double qt = sqrt(-t * (-t / cnts->t_min + 1.));

	TComplex sum = 0.;
	for (int j = 1; j <= cgcMaxOrder; j++)
	{
		double F = F_cal(j, qt, lambda, m0sq);
		sum += F*F * T_lxg_n(j, t);
	}

	return Abs_fac * sum;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------- FULL AMPLITUDE ----------------------------------------

TComplex IslamModel::Amp(double t) const
{	
#ifdef DEBUG
	printf(">> IslamModel::amp, mode = %i\n", mode);
#endif

	switch (mode)
	{
		case mDiff:		return T_diff(t);
		case mCore:		return T_core(t);
		case mQuark:	return T_quark(t);
		case mDiffCore:	return T_diff(t) + T_core(t);
		case mFull:	 	return T_diff(t) + T_core(t) + T_quark(t);
		default:
			printf("ERROR in IslamModel::Amp > unknown mode %i\n", mode);
			return 0.;
	}
}


//----------------------------------------------------------------------------------------------------
//---------------------------------------- PROFILE FUNCTIONS --------------------------------------

TComplex IslamModel::Amp_J0(double t, double *par, const void *vobj)
{
	const IslamModel *obj = (IslamModel *) vobj;
	const double &b = par[0];	

	return obj->Amp(t) * TMath::BesselJ0(b*sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::Prf(double b_fm) const
{
	double b = b_fm / cnts->hbarc;	// b in GeV^-1
	double par[] = { b };

	TComplex I = ComplexIntegrate(Amp_J0, par, this, upper_bound_t, 0., 0., precision_t,
		integ_workspace_size_t, integ_workspace_t, "IslamModel::Prf");
	return I / 4. / cnts->p_cms / cnts->sqrt_s;
}
