// HEADER

#include "interface/PPPModel.h"
#include "interface/Math.h"
#include "interface/Constants.h"

using namespace std;
using namespace Elegent;

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

void PPPModel::Init(unsigned char m)
{
	// physics parameters
	s0 = 1;	// in GeV^2
	
	tag = "";
	//						Delta		c			a'		r^2
	//						1			1			GeV^2	GeV^2
	if (m == m2P)
	{
		// two pomerons
		mode = m2P;
		SetTrajectory(pom1, 0.08590,	53.18,		0.360,	9.595,	 s0);
		SetTrajectory(pom2, 0.14437,	6.87,		0.082,	4.765,	 s0);
		SetTrajectory(oder, -0.2707,	1.8134,		0.029,	1.159,	s0);
		SetTrajectory(regf, -0.3100,	188.51,		0.84,	41.424, s0);
		SetTrajectory(rego, -0.5300,	-171.36,	0.93,	2.621,	s0);
		tag = "ppp2";
		name = "PPP2";
		label = "Petrov-Predazzi-Prokudin, 2 pomerons";
	}
	
	if (m == m3P)
	{
		// three pomerons
		mode = m3P;
		SetTrajectory(pom1, 0.0578,		53.007,	 0.5596, 6.3096, s0);
		SetTrajectory(pom2, 0.1669,		9.6762,		0.2733,	3.1097, s0);
		SetTrajectory(pom3, 0.2032,		1.6654,		0.0937,	2.4771, s0);
		SetTrajectory(oder, 0.1920,		0.0166,		0.048,	0.1398, s0);
		SetTrajectory(regf, -0.31,		191.69,		0.84,	31.593, s0);
		SetTrajectory(rego, -0.53,		-174.18,	0.93,	7.467,	 s0);
		tag = "ppp3";
		name = "PPP3";
		label = "Petrov-Predazzi-Prokudin, 3 pomerons";
	}

	if (tag.empty())
		printf("ERROR in PPPModel::Init > Unknown mode %u.\n", m);
	
	// integration parameters
	precision = 1E-14;
	upper_bound = 60;
}

//----------------------------------------------------------------------------------------------------

void PPPModel::Print() const
{
	printf(">> PPPModel::Print\n");
	printf("\tmode: %s\n", GetModeString().c_str());

	printf("\tpom1: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", pom1.D, pom1.c, pom1.ap, pom1.r2);
	printf("\tpom2: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", pom2.D, pom2.c, pom2.ap, pom2.r2);
	if (mode == m3P)
		printf("\tpom3: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", pom3.D, pom3.c, pom3.ap, pom3.r2);
	printf("\toder: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", oder.D, oder.c, oder.ap, oder.r2);
	printf("\tregf: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", regf.D, regf.c, regf.ap, regf.r2);
	printf("\trego: delta=%.4f, c=%.4f, a'=%.4f, r^2=%.4f\n", rego.D, rego.c, rego.ap, rego.r2);
}

//----------------------------------------------------------------------------------------------------

string PPPModel::GetModeString() const
{
	switch (mode)
	{
		case m2P: return "2 pomerons"; 
		case m3P: return "3 pomerons"; 
		default: return "unknown";
	}
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::tr_eik(Trajectory t, double , double b) const
{
	// delta+-(s, b) ... eq (11), (12) without the leading i
	// the s-dependence is now computed during initialization
	// therefore any cnts->s change after init() won't have any affect
	return t.gamma * exp(- b * b / t.rho2) / (4. * cnts->pi * t.rho2);
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::prf0(double b) const
{
	TComplex delta = i*tr_eik(pom1, cnts->s, b) + i*tr_eik(pom2, cnts->s, b) + i*tr_eik(regf, cnts->s, b);

	if (cnts->pMode == cnts->mPP)
		delta += tr_eik(oder, cnts->s, b) + tr_eik(rego, cnts->s, b);
	else
		delta -= tr_eik(oder, cnts->s, b) + tr_eik(rego, cnts->s, b);

	if (mode == m3P)
		delta += i*tr_eik(pom3, cnts->s, b);

	return (TComplex::Exp(2.*i*delta) - 1.) / 2. / i;
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::Prf(double b) const
{
	return prf0(b / cnts->hbarc); 
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::prf_J0(double *b, double *t, const void *obj)
{
	return ((PPPModel *)obj)->prf0(b[0]) * b[0] * TMath::BesselJ0(b[0] * sqrt(-t[0]));
}

//----------------------------------------------------------------------------------------------------

TComplex PPPModel::Amp(double t) const
{
	return 2.*cnts->p_cms*cnts->sqrt_s * CmplxInt(this, prf_J0, 0, upper_bound, &t, precision);
}