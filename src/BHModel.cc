/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#include "interface/BHModel.h"
#include "interface/Math.h"
#include "interface/Constants.h"

using namespace Elegent;

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

	// TODO: which is correct?
	// [1] : Ng = (6.-epsilon)*(5.-epsilon)*(4.-epsilon)*(3.-epsilon)*(2.-epsilon)*(1.-epsilon) /5./4./3./2. / 2.;
	// [2] : Ng = 3./2. * (5.-epsilon)*(4.-epsilon)*(3.-epsilon)*(2.-epsilon)*(1.-epsilon) /5./4./3./2.;
	Ng = (6.-epsilon)*(5.-epsilon)*(4.-epsilon)*(3.-epsilon)*(2.-epsilon)*(1.-epsilon) /5./4./3./2. / 2.;
		
	a0 = -41.1;
	a1 = -487.5;
	a2 = -600.;
	a3 = 600.;
	a4 = 487.5;
	a5 = 41.1;

	b0 = -9.;
	b1 = -225.;
	b2 = -900.;
	b3 = -900.;
	b4 = -225.;
	b5 = -9.;

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

	// integration settings
	upper_bound = 50.;
	precision = 1E-12;

	// precompute mu_qg
	mu_qg = sqrt(mu_qq*mu_gg);
		
	// precompute sigma_gg, Eq. (B5) in [1]: below the Cp_gg stands for C'_gg
	sigma_gg = Cp_gg * Sigma_gg * Ng*Ng * TComplex(sumR(cnts->s), sumI(cnts->s));

	// precompute sigma_qq, Eq. (B9) in [1]
	sigma_qq = Sigma_gg * (C + C_even_regge * m0/sqrt(cnts->s) * TComplex::Exp(i * cnts->pi / 4.));

	// precompute sigma_qg, Eq. (B10) in [1]
	sigma_qg = Sigma_gg * C_qg_log * TComplex(log(cnts->s/s0), -cnts->pi/2.);
	
	// precompute sigma_odd, Eq. (B12) in [1] - the factor in front of W(b, mu_odd)
	// plus additional factor (-i) to match the normalization used in chi_without_i
	sigma_odd = -i * C_odd * Sigma_gg * m0 / cnts->sqrt_s * TComplex::Exp(i * cnts->pi / 4.);
}

//----------------------------------------------------------------------------------------------------

void BHModel::Print() const
{
	printf(">> BHModel::Print\n");

	printf("\tcommon\n");
	printf("\t\ts0 = %E\n", s0);
	printf("\t\tm0 = %E\n", m0);
	printf("\t\tSigma_gg = %E\n", Sigma_gg);

	printf("\tsigma_gg\n");
	printf("\t\tCp_gg = %E\n", Cp_gg);
	printf("\t\tNg = %E\n", Ng);
	printf("\t\tepsilon = %E\n", epsilon);
	printf("\t\ta0 = %E, a1 = %E, a2 = %E, a3 = %E, a4 = %E, a5 = %E\n", a0, a1, a2, a3, a4, a5);
	printf("\t\tb0 = %E, b1 = %E, b2 = %E, b3 = %E, b4 = %E, b5 = %E\n", b0, b1, b2, b3, b4, b5);
	printf("\t\tsigma_gg: Re = %E, Im = %E\n", sigma_gg.Re(), sigma_gg.Im());

	printf("\tsigma_qg\n");
	printf("\t\tC_qg_log = %E\n", C_qg_log);
	printf("\t\tsigma_qg: Re = %E, Im = %E\n", sigma_qg.Re(), sigma_qg.Im());

	printf("\tsigma_qq\n");
	printf("\t\tC = %E\n", C);
	printf("\t\tC_even_regge = %E\n", C_even_regge);
	printf("\t\tsigma_qq: Re = %E, Im = %E\n", sigma_qq.Re(), sigma_qq.Im());

	printf("\tsigma_odd\n");
	printf("\t\tC_odd = %E\n", C_odd);
	printf("\t\tsigma_odd: Re = %E, Im = %E\n", sigma_odd.Re(), sigma_odd.Im());
	
	printf("\tmu's\n");
	printf("\t\tmu_gg = %E\n", mu_gg);
	printf("\t\tmu_qq = %E\n", mu_qq);
	printf("\t\tmu_qg = %E\n", mu_qg);
	printf("\t\tmu_odd = %E\n", mu_odd);
}

//----------------------------------------------------------------------------------------------------

double BHModel::sumR1(double s) const
{
	return
		(a0-b0/(-epsilon))/(-epsilon)+(a1-b1/(1-epsilon))/(1-epsilon)
		+(a2-b2/(2-epsilon))/(2-epsilon)+(a3-b3/(3-epsilon))/(3-epsilon)
		+(a4-b4/(4-epsilon))/(4-epsilon)+(a5-b5/(5-epsilon))/(5-epsilon)
		-(a0-b0/(-epsilon))/(-epsilon)*exp(epsilon*log(s/m0/m0))*cos(epsilon*cnts->pi/2)
		-(a1-b1/(1-epsilon))/(1-epsilon)*exp((epsilon-1)*log(s/m0/m0))*cos((epsilon-1)*cnts->pi/2)
		-(a2-b2/(2-epsilon))/(2-epsilon)*exp((epsilon-2)*log(s/m0/m0))*cos((epsilon-2)*cnts->pi/2)
		-(a3-b3/(3-epsilon))/(3-epsilon)*exp((epsilon-3)*log(s/m0/m0))*cos((epsilon-3)*cnts->pi/2)
		-(a4-b4/(4-epsilon))/(4-epsilon)*exp((epsilon-4)*log(s/m0/m0))*cos((epsilon-4)*cnts->pi/2)
		-(a5-b5/(5-epsilon))/(5-epsilon)*exp((epsilon-5)*log(s/m0/m0))*cos((epsilon-5)*cnts->pi/2);
}

//----------------------------------------------------------------------------------------------------

double BHModel::sumR2(double s) const
{
	return
		-(b0/(-epsilon))*exp(epsilon*log(s/m0/m0))*cos(epsilon*cnts->pi/2)*log(m0*m0/s)
		-(b1/(1-epsilon))*exp((epsilon-1)*log(s/m0/m0))*cos((epsilon-1)*cnts->pi/2)*log(m0*m0/s)
		-(b2/(2-epsilon))*exp((epsilon-2)*log(s/m0/m0))*cos((epsilon-2)*cnts->pi/2)*log(m0*m0/s)
		-(b3/(3-epsilon))*exp((epsilon-3)*log(s/m0/m0))*cos((epsilon-3)*cnts->pi/2)*log(m0*m0/s)
		-(b4/(4-epsilon))*exp((epsilon-4)*log(s/m0/m0))*cos((epsilon-4)*cnts->pi/2)*log(m0*m0/s)
		-(b5/(5-epsilon))*exp((epsilon-5)*log(s/m0/m0))*cos((epsilon-5)*cnts->pi/2)*log(m0*m0/s)
		-(b0/(-epsilon))*exp(epsilon*log(s/m0/m0))*sin(epsilon*cnts->pi/2)*cnts->pi/2
		-(b1/(1-epsilon))*exp((epsilon-1)*log(s/m0/m0))*sin((epsilon-1)*cnts->pi/2)*cnts->pi/2
		-(b2/(2-epsilon))*exp((epsilon-2)*log(s/m0/m0))*sin((epsilon-2)*cnts->pi/2)*cnts->pi/2
		-(b3/(3-epsilon))*exp((epsilon-3)*log(s/m0/m0))*sin((epsilon-3)*cnts->pi/2)*cnts->pi/2
		-(b4/(4-epsilon))*exp((epsilon-4)*log(s/m0/m0))*sin((epsilon-4)*cnts->pi/2)*cnts->pi/2
		-(b5/(5-epsilon))*exp((epsilon-5)*log(s/m0/m0))*sin((epsilon-5)*cnts->pi/2)*cnts->pi/2;
}

//----------------------------------------------------------------------------------------------------

double BHModel::sumR(double s) const
{
	return sumR1(s) + sumR2(s);
}

//----------------------------------------------------------------------------------------------------

double BHModel::sumI1(double s) const
{
	return
		(a0-b0/(-epsilon))/(-epsilon)*exp(epsilon*log(s/m0/m0))*sin(epsilon*cnts->pi/2)
		+(a1-b1/(1-epsilon))/(1-epsilon)*exp((epsilon-1)*log(s/m0/m0))*sin((epsilon-1)*cnts->pi/2)
		+(a2-b2/(2-epsilon))/(2-epsilon)*exp((epsilon-2)*log(s/m0/m0))*sin((epsilon-2)*cnts->pi/2)
		+(a3-b3/(3-epsilon))/(3-epsilon)*exp((epsilon-3)*log(s/m0/m0))*sin((epsilon-3)*cnts->pi/2)
		+(a4-b4/(4-epsilon))/(4-epsilon)*exp((epsilon-4)*log(s/m0/m0))*sin((epsilon-4)*cnts->pi/2)
		+(a5-b5/(5-epsilon))/(5-epsilon)*exp((epsilon-5)*log(s/m0/m0))*sin((epsilon-5)*cnts->pi/2)
		-(b0/(-epsilon))*exp(epsilon*log(s/m0/m0))*cos(epsilon*cnts->pi/2)*cnts->pi/2
		-(b1/(1-epsilon))*exp((epsilon-1)*log(s/m0/m0))*cos((epsilon-1)*cnts->pi/2)*cnts->pi/2
		-(b2/(2-epsilon))*exp((epsilon-2)*log(s/m0/m0))*cos((epsilon-2)*cnts->pi/2)*cnts->pi/2
		-(b3/(3-epsilon))*exp((epsilon-3)*log(s/m0/m0))*cos((epsilon-3)*cnts->pi/2)*cnts->pi/2
		-(b4/(4-epsilon))*exp((epsilon-4)*log(s/m0/m0))*cos((epsilon-4)*cnts->pi/2)*cnts->pi/2
		-(b5/(5-epsilon))*exp((epsilon-5)*log(s/m0/m0))*cos((epsilon-5)*cnts->pi/2)*cnts->pi/2;
}

//----------------------------------------------------------------------------------------------------

double BHModel::sumI2(double s) const
{
	return
		(b0/(-epsilon))*exp(epsilon*log(s/m0/m0))*sin(epsilon*cnts->pi/2)*log(m0*m0/s)
		+(b1/(1-epsilon))*exp((epsilon-1)*log(s/m0/m0))*sin((epsilon-1)*cnts->pi/2)*log(m0*m0/s)
		+(b2/(2-epsilon))*exp((epsilon-2)*log(s/m0/m0))*sin((epsilon-2)*cnts->pi/2)*log(m0*m0/s)
		+(b3/(3-epsilon))*exp((epsilon-3)*log(s/m0/m0))*sin((epsilon-3)*cnts->pi/2)*log(m0*m0/s)
		+(b4/(4-epsilon))*exp((epsilon-4)*log(s/m0/m0))*sin((epsilon-4)*cnts->pi/2)*log(m0*m0/s)
		+(b5/(5-epsilon))*exp((epsilon-5)*log(s/m0/m0))*sin((epsilon-5)*cnts->pi/2)*log(m0*m0/s);
}

//----------------------------------------------------------------------------------------------------

double BHModel::sumI(double s) const
{
	return sumI1(s) + sumI2(s);
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
	// TODO: why the one half?

	// Eqs. (B1) without the leading i factor and Eq. (B12) in [1]
	return (
			sigma_gg * W(b, mu_gg) + sigma_qg * W(b, mu_qg) + sigma_qq * W(b, mu_qq)
			+ ((cnts->pMode == cnts->mAPP) ? +1. : -1.) * sigma_odd * W(b, mu_odd)
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

TComplex BHModel::prf0_J0(double *b, double *q, const void *obj)
{
	return b[0] * TMath::BesselJ0(b[0]*q[0]) * ((BHModel *)obj)->prf0(b[0]);
} 

//----------------------------------------------------------------------------------------------------

TComplex BHModel::Amp(double t) const
{
	double q = sqrt(-t);

	// from Eqs. (A11) and (A12)
	return	i * cnts->p_cms * cnts->sqrt_s * CmplxInt(this, prf0_J0, 0., upper_bound, &q, precision);
}
