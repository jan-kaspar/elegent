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

#include "interface/Constants.h"
#include "interface/ModelFactory.h"
#include "interface/CoulombInterference.h"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"

using namespace Elegent;
using namespace std;

//----------------------------------------------------------------------------------------------------

/// hadronic model for fitting
class HadronicFitModel : public Model
{
  public:
	/// modulus parameters (low |t|)
	double a, b1, b2, b3;

	/// parameters for blending the low-|t| (variable) and high-|t| (fixed) modulus
	/// the interval (t1, t2) corresponds to (-3 sigma, +3 sigma)
	double t1, t2;

	/// phase parameters
	double p0;

    HadronicFitModel() :
		a(0.), b1(0.), b2(0.), b3(0.),
		t1(0.1), t2(0.4),
		p0(0.)
	{
	}

	/// elastic amplitude
    virtual TComplex Amp(double t) const
	{
		// low-|t| modulus
		double bPol = 0., tPow = t;
		bPol += b1 * tPow; tPow *= t;
		bPol += b2 * tPow; tPow *= t;
		bPol += b3 * tPow; tPow *= t;
		const double m1 = a * exp(bPol);

		// high-|t| modulus
		const double P0 = +6.027359E+02;
		const double P1 = +9.732869E+01;
		const double P2 = -2.031860E+01;
		const double P3 = +5.145647E+00;
		const double P4 = -9.061766E+00;
		const double P5 = -1.245398E+01;
		const double P6 = +2.673420E+01;
		const double P7 = -6.135371E+00;
		const double x = -t;
		const double dsdt2 = (P0 + P1*x) * exp(P2*x + P3*x*x + P4*x*x*x) + (P5 + P6*x) * exp(P7*x);
		const double m2 = sqrt(dsdt2 / cnts->sig_fac);

		// full (blended) modulus
		const double t_avg = (t2 + t1)/2., t_si = (t2 - t_avg) / 3.;
		const double f = (TMath::Erf( (-t - t_avg) / t_si / sqrt(2.) ) + 1.) / 2.;
		const double m = m1*(1.-f) + m2*f;

		// phase
		const double p = p0;

		// full amplitude
		return m * TComplex::Exp(i * p);
	}

	virtual void Init() {}

	virtual void Print() const
	{
		printf("hfm: a = %.3E, b1 = %.3E, b2 = %.3E, b3 = %.3E, p0 = %.3E\n", a, b1, b2, b3, p0);
	}

	virtual TComplex Prf(double) const { return 0; }
};

//----------------------------------------------------------------------------------------------------

/// instance of the model object used in the fit
HadronicFitModel *hfm;

/// fit function
double f_fit_imp(double x[], double par[])
{
	// transfer parameters to hadronic fit model
	hfm->a = par[0] * 1E9;
	hfm->b1 = par[1];
	hfm->b2 = par[2];
	hfm->b3 = par[3];
	hfm->p0 = par[4];

	// evaluate differential cross-section
	const double t = -x[0];
	return cnts->sig_fac * coulomb->Amp(t).Rho2();
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// initialize calculation engine: for proton-proton scattering at sqrt(s) = 13 TeV
	Constants::Init(2*6500., cnts->mPP);
    cnts->Print();

	// initialize CNI calculation
	coulomb->mode = CoulombInterference::mKL;	// formula by Kundrat-Lokajicek
	coulomb->ffType = coulomb->ffPuckett;
	coulomb->precision = 1E-3;
	coulomb->Print();

	// switch on optimisation in CNI calculation
	coulomb->InitIFunctionInterpolator(20., 300);
	coulomb->useIFunctionInterpolator = true;

	// prepare output file
	TFile *f_out = TFile::Open("example2.root", "recreate");

	// build sample data - differential cross-section
	//  - for example model by Petrov et al. (with 3 Pomerons)
	//  - for example 1% uncertainty
	PPPModel *ppp3 = new PPPModel();
	ppp3->Configure(PPPModel::v3P);
	ppp3->Init();
	model = ppp3;

	printf("* rho in simulation: %.3f\n", model->Amp(0).Re() / model->Amp(0).Im());

	printf("* running simulation\n");
	TGraphErrors *g_dsdt = new TGraphErrors();
	g_dsdt->SetName("g_dsdt");

	for (double mt = 1E-4; mt < 0.1; mt *= 1.2)
	{
		const double t = -mt;	// four-momentum transfer squared, in GeV^2
		const TComplex A = coulomb->Amp(t);		// scattering amplitude
		const double dsdt = cnts->sig_fac * A.Rho2();	// differential cross-section

		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, mt, dsdt);
		g_dsdt->SetPointError(idx, 0., dsdt * 0.01);
	}

	g_dsdt->Write();

	// initialise fit model
	hfm = new HadronicFitModel();
	model = hfm;

	// initialise fit function
	TF1 *f_fit = new TF1("f_fit", f_fit_imp, 1E-5, 1., 5);
	f_fit->SetNpx(1000);

	// start point for minimisation
	f_fit->SetParameter(0, 1.857);					// a / 1E9
	f_fit->SetParameter(1, 10.1);					// b1
	f_fit->SetParameter(2, 0.);						// b2
	f_fit->SetParameter(3, 0.);						// b3
	f_fit->SetParameter(4, M_PI/2. - atan(0.14));	// rho = 0.14

	// run minimisation
	printf("* running fit\n");
	g_dsdt->Fit(f_fit);

	// print result
	printf("* rho in fit: %.3f\n", 1./tan(f_fit->GetParameter(4)));

	// save fit graph
	TGraph *g_dsdt_fit = new TGraph();
	g_dsdt_fit->SetName("g_dsdt_fit");
	g_dsdt_fit->SetLineColor(2);
	for (double mt = 1E-4; mt < 0.1; mt *= 1.2)
	{
		const double t = -mt;	// four-momentum transfer squared, in GeV^2
		const TComplex A = coulomb->Amp(t);		// scattering amplitude
		const double dsdt = cnts->sig_fac * A.Rho2();	// differential cross-section

		int idx = g_dsdt_fit->GetN();
		g_dsdt_fit->SetPoint(idx, mt, dsdt);
	}
	g_dsdt_fit->Write();

	// clean up
	delete f_out;

	return 0;
}
