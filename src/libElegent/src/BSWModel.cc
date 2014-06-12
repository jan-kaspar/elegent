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

#include "interface/BSWModel.h"
#include "interface/Constants.h"

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

BSWModel::BSWModel()
{
	fullLabel.name = "Bourrely et al."; shortLabel.name = "bourrely";

	integ_workspace_initialized = false;
}

//----------------------------------------------------------------------------------------------------

BSWModel::~BSWModel()
{
	if (integ_workspace_initialized)
	{
		gsl_integration_workspace_free(integ_workspace_b);
		gsl_integration_workspace_free(integ_workspace_t);
	}
}

//----------------------------------------------------------------------------------------------------

void BSWModel::Configure(BSWModel::ModeType _mode, bool _presampled)
{
	mode = _mode;
	presampled = _presampled;

	// set labels
	fullLabel.variant = ""; shortLabel.variant = "";
	fullLabel.version = "Eur. Phys. J. C28 (2003) 97-105"; shortLabel.version = "03";
	if (mode == mPomReg)
		{ fullLabel.mode = "full"; shortLabel.mode = "full"; }
	if (mode == mPom)
		{ fullLabel.mode = "Pomeron"; shortLabel.mode = "pom"; }
	if (mode == mReg)
		{ fullLabel.mode = "Reggeon"; shortLabel.mode = "reg"; }

	// ambiguity constants
	k_u = 0; k_lnu = 0;	// only this combination reproduces high-energy behaviour of sigma_tot(s) and rho(s)
	regge_fac(0., 1.);	// only way to reproduce results at low energies
}

//----------------------------------------------------------------------------------------------------

void BSWModel::Init()
{
	// Eq. (6) in [3]
	c = 0.167;
	cp = 0.748;

	// Table 1 in [3]
	m1 = 0.577; m1sq = m1*m1;	// 0.333
	m2 = 1.719; m2sq = m2*m2;	// 2.955
	f = 6.971;
	a = 1.858;	asq = a*a;		// 3.452
	 
	// Table 2 in [3], b's given in text of section 3 in [3]
	A2.Init(	-24.269,	0.,		0.357,	1.,		+1);
	omega.Init( -167.329,	0.,		0.323,	0.795,	-1);
	rho.Init(	124.919,	8.54,	0.320,	1.,		-1);

	// signs according to interaction mode
	if (cnts->pMode == cnts->mPP)
	{
		A2.sign = +1.; omega.sign = rho.sign = +1.;
	} else {
		A2.sign = +1.; omega.sign = rho.sign = -1.;
	}

	// integration parameters
	upper_bound_t = -500.; precision_t = 1E-3;
	upper_bound_b = 50.; precision_b = 1E-1;

	// prepare integration workspace
	if (!integ_workspace_initialized)
	{
		integ_workspace_size_b = 100;
		integ_workspace_b = gsl_integration_workspace_alloc(integ_workspace_size_b);

		integ_workspace_size_t = 100;
		integ_workspace_t = gsl_integration_workspace_alloc(integ_workspace_size_b);

		integ_workspace_initialized = true;
	}

	// save value of S0(0)
	S00 = S0(0.);
	if (presampled)
		BuildSample(10001);
}

//----------------------------------------------------------------------------------------------------

void BSWModel::Print() const
{
	printf(">> BSWModel::Print\n");
	printf("\t%s\n", CompileFullLabel().c_str());
	printf("\tmode = %u\n", mode);
	printf("\tc=%.3f, c'=%.3f, m1=%.3f, m2=%.3f, f=%.3f, a=%.3f\n", c, cp, m1, m2, f, a);
	printf("\tA2   : C=%.3f, b=%.3f, alpha=%.3f, aplha'=%.3f, signature=%+.0f, sign=%+0.f\n",
		A2.C, A2.b, A2.a, A2.ap, A2.signature, A2.sign);
	printf("\tomega: C=%.3f, b=%.3f, alpha=%.3f, aplha'=%.3f, signature=%+.0f, sign=%+0.f\n",
		omega.C, omega.b, omega.a, omega.ap, omega.signature, omega.sign);
	printf("\trho  : C=%.3f, b=%.3f, alpha=%.3f, aplha'=%.3f, signature=%+.0f, sign=%+0.f\n",
		rho.C, rho.b, rho.a, rho.ap, rho.signature, rho.sign);
	printf("\tregge_fac = %+.0f %+.0fi\n", regge_fac.Re(), regge_fac.Im());
	printf("\tk_u = %i, k_lnu = %i\n", k_u, k_lnu);

	printf("\n");

	printf("\tpresampled = %u\n", presampled);
	if (presampled)
		printf("\t\tsample size = %u, db = %.1E\n", data_N, data_db);
	
	printf("\n");

	printf("\tintegration parameters:\n");
	printf("\t\tt: upper bound = %.1E, precision = %.1E\n", upper_bound_t, precision_t);
	printf("\t\tb: upper bound = %.1E, precision = %.1E\n", upper_bound_b, precision_b);
}

//----------------------------------------------------------------------------------------------------

double BSWModel::Ft(double t) const
{
	/// Eq. (4) in [3]
	double G = 1. / (1. - t/m1sq) / (1. - t/m2sq);
	return f*G*G*	(asq + t) / (asq - t);
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::Rt(Trajectory tr, double t) const
{
	/// Eq. (7) in [3]
	/// NOTE: the equation is missing signs that make the difference between pp and app, corrected
	/// version below.

	/// NOTE: as s0 = 1 GeV^2, it is simply left out of the formula below.

	double alpha = tr.a + tr.ap * t;
	return tr.sign * tr.C * exp(tr.b * t + cnts->ln_s * alpha) *
		(1. + tr.signature * TComplex::Exp(-i*cnts->pi*alpha));
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::R0t(double t) const
{
	/// Sum over all allowed Regge trajectories, as defined in the text below Eq. (7) in [3].
	/// version below.

	return Rt(A2, t) + Rt(omega, t) + Rt(rho, t);
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::S0(double t) const
{
	/// Eq. (3) in [3]

	// the s term
	TComplex term_s = pow(cnts->s, c) / pow(cnts->ln_s, cp);

	// the u term
	double u = 4.*cnts->proton_mass*cnts->proton_mass - cnts->s - t;
	
	// u = -|u| exp( i * (2 k_u pi - pi) ) ==> ln u = ln |u| + i * (2 k_u pi - pi)
	TComplex Lnu = TComplex(log(fabs(u)), (2.* k_u - 1.) * cnts->pi);

	double Lnu_rho2 = Lnu.Rho2();
	double Lnu_theta = atan2(Lnu.Im(), Lnu.Re());	 // atan2 results in (-pi, +pi)

	// ln u = sqrt(Lu_rho2) exp(i * (Lnu_theta + 2 k_lnu pi))
	TComplex LnLnu = TComplex(0.5*log(Lnu_rho2), Lnu_theta + 2. * k_lnu * cnts->pi);

	TComplex term_u = TComplex::Exp(c * Lnu) / TComplex::Exp(cp * LnLnu);
	
#ifdef DEBUG
	printf(">> BSWModel::S0\n");
	printf("\ts=%E, u=%E\n", cnts->s, u);
	printf("\tLn u=%E +i%E\n", Lnu.Re(), Lnu.Im());
	printf("\tLn Ln u=%E +i%E\n", LnLnu.Re(), LnLnu.Im());
	printf("\tterm_s: %E + i%E\n", term_s.Re(), term_s.Im());
	printf("\tterm_u: %E + i%E\n", term_u.Re(), term_u.Im());
#endif

	return term_s + term_u;
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::Omega0t(double t) const
{
	/// Eq. (2) in [3] written in t-space (cf. Eq. (7) in [1])
	/// S00 = S0(0) instead of S0(t) is used here, valid for high s only !!

	switch (mode)
	{
		case mPomReg:
			return S00*Ft(t) + R0t(t) / cnts->s / regge_fac;
		case mPom:
			return S00*Ft(t);
		case mReg:
			return R0t(t);
	}
	
	return TComplex(0, 0);
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::Omega0t_J0(double t, double *par, const void *vobj)
{
	const BSWModel *obj = (BSWModel *) vobj;
	const double &b = par[0];

	return obj->Omega0t(t) * TMath::BesselJ0(b * sqrt(-t));
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::Omega0b(double b) const
{
	// the 1/2 factor is consequence of dt integration (instead of q dq)
	double par[] = { b };
	return 0.5 * ComplexIntegrate(Omega0t_J0, par, this, upper_bound_t, 0., precision_t,
		integ_workspace_size_t, integ_workspace_t, "BSWModel::Omega0b");
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::prf0(double b) const
{
	return 1. - TComplex::Exp( - Omega0b(b) );
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::Prf(double b) const
{
	return prf0(b / cnts->hbarc) * i / 2.;
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::prf0_J0(double b, double *par, const void *vobj)
{
	const BSWModel *obj = (BSWModel *) vobj;
	const double &q = par[0];

	TComplex prf0_v = (obj->presampled) ? obj->SampleEval(b) : obj->prf0(b);
	return prf0_v  * b * TMath::BesselJ0(b * q);
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::Amp(double t) const
{
	double q = sqrt(-t);
	double par[] = { q };

	TComplex I = ComplexIntegrate(prf0_J0, par, this, 0., upper_bound_b, precision_b,
		integ_workspace_size_b, integ_workspace_b, "BSWModel::Amp");
	
	return i * cnts->p_cms * cnts->sqrt_s * I;
}

//----------------------------------------------------------------------------------------------------

void BSWModel::BuildSample(unsigned int samples)
{
	printf(">> BSWModel::BuildSample > Building %u samples...\n", samples);

	data_re.clear();
	data_re.reserve(samples);
	data_im.clear();
	data_im.reserve(samples);
#ifdef DEBUG
	data_b.clear();
	data_b.reserve(samples);
#endif

	double db = upper_bound_b / (samples - 1);
	data_db = db;
	data_N = samples;

	double b = 0.;
	for (unsigned int i = 0; i < samples; i++, b += db)
	{
		TComplex v = prf0(b);

		//printf("v=%.5f: re=%E, im=%E\n", b, v.Re(), v.Im());
#ifdef DEBUG
		data_b.push_back(b);
#endif

		data_re.push_back(v.Re());
		data_im.push_back(v.Im());
	}
}

//----------------------------------------------------------------------------------------------------

TComplex BSWModel::SampleEval(double b) const
{
	unsigned int idx = (int)(b / data_db);
	
	if (idx + 1 > data_N - 1)
		return TComplex(0, 0);

	double f = b/data_db - idx;

	return TComplex(
		(data_re[idx+1] - data_re[idx])*f + data_re[idx],
		(data_im[idx+1] - data_im[idx])*f + data_im[idx]
	);
}
