#include "interface/Constants.h"
#include "interface/CoulombInterference.h"

#include "HadronicFitModel.h"

#include "TFile.h"
#include "TGraph.h"

#include <string>
#include <chrono>
#include <iostream>

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

int main()
{
	// init Elegent
	Constants::Init(2*6500, cnts->mPP);
    cnts->Print();

	coulomb->mode = CoulombInterference::mKL;
	coulomb->ffType = coulomb->ffPuckett;
	coulomb->precision = 1E-5;
	coulomb->Print();

	// init hadronic model
	HadronicFitModel *hfm = new HadronicFitModel();

	hfm->a = 1.84E9;
	hfm->b1 = 10.2;
	hfm->p0 = M_PI/2. - atan(0.12);

	hfm->t1 = 0.2;
	hfm->t2 = 0.5;

	hfm->modulusHighTVariant = 2;
	hfm->phaseMode = HadronicFitModel::pmConstant;

	model = hfm;

	// prepare output
	TFile *f_out = TFile::Open("coulomb_I_function_interpolation_t.root", "recreate");

	// timing variables
	chrono::time_point<chrono::system_clock> t_start, t_end;
	chrono::duration<double> t_duration;

	// prepare interpolator
	cout << "* initiating interpolator" << endl;
	t_start = chrono::system_clock::now();
	coulomb->InitIFunctionInterpolator(20., 300);
	t_end = chrono::system_clock::now();

	t_duration = t_end - t_start;
	cout << "    done in " << t_duration.count() << " s" << endl;

	// define t sampling
	double mt_min = 1E-5, mt_max = 0.2, mt_step = 2E-4;

	// evaluate t-distribution, I(t, t'): numerical integration
	coulomb->useIFunctionInterpolator = false;
	TGraph *g_dsdt_integration = new TGraph();
	cout << "* evaluation t-distribution, I: numerical integration" << endl;
	t_start = chrono::system_clock::now();
	for (double mt = mt_min; mt < mt_max; mt += mt_step)
	{
		double dsdt = cnts->sig_fac * coulomb->Amp(-mt).Rho2();

		int idx = g_dsdt_integration->GetN();
		g_dsdt_integration->SetPoint(idx, mt, dsdt);
	}
	g_dsdt_integration->Write("g_dsdt_integration");
	t_end = chrono::system_clock::now();
	t_duration = t_end - t_start;
	cout << "    done in " << t_duration.count() << " s" << endl;

	// evaluate t-distribution, I(t, t'): interpolation
	coulomb->useIFunctionInterpolator = true;
	TGraph *g_dsdt_interpolation = new TGraph();
	cout << "* evaluation t-distribution, I: interpolation" << endl;
	t_start = chrono::system_clock::now();
	for (double mt = mt_min; mt < mt_max; mt += mt_step)
	{
		double dsdt = cnts->sig_fac * coulomb->Amp(-mt).Rho2();

		int idx = g_dsdt_interpolation->GetN();
		g_dsdt_interpolation->SetPoint(idx, mt, dsdt);
	}
	g_dsdt_interpolation->Write("g_dsdt_interpolation");
	t_end = chrono::system_clock::now();
	t_duration = t_end - t_start;
	cout << "    done in " << t_duration.count() << " s" << endl;

	// clean up
	delete f_out;

	return 0;
}
