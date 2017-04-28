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
#include "interface/CoulombInterference.h"
#include "interface/ModelFactory.h"

#include "TFile.h"
#include "TGraph.h"

using namespace Elegent;
using namespace std;

//----------------------------------------------------------------------------------------------------

void TestIslam()
{
	vector<double> energies;
	vector<Constants::ParticleMode> pModes;
	energies.push_back(1960); pModes.push_back(Constants::mAPP);
	energies.push_back(7000); pModes.push_back(Constants::mPP);
	energies.push_back(8000); pModes.push_back(Constants::mPP);
	energies.push_back(13000); pModes.push_back(Constants::mPP);

	//vector<int> modes = { 0, 1, 2, 3, 4 };

	for (unsigned int ei = 0; ei < energies.size(); ei++)
	{
		char buf[100];
		sprintf(buf, "debugger_%.0f.root", energies[ei]);
		TFile *f_out = new TFile(buf, "recreate");
	
		Constants::Init(energies[ei], pModes[ei]);

		printf(">> s = %.3E\n", cnts->s);

		for (unsigned int mode = 0; mode <= 4; mode++)
		{
			IslamModel2015 *im = new IslamModel2015;
			im->Configure((IslamModel2015::ModeType) mode);
			im->Init();

			TGraph *g_dsdt = new TGraph();
			sprintf(buf, "g_dsdt_%i", mode);
			g_dsdt->SetName(buf);
			g_dsdt->SetTitle(im->fullLabel.mode.c_str());

			for (double mt = 1e-4; mt < 2.5; mt += 0.001)
			{
				double dsdt = cnts->sig_fac * im->Amp(-mt).Rho2();
				//printf("t = %E: dsdt = %E\n", mt, dsdt);
		
				int idx = g_dsdt->GetN();
				g_dsdt->SetPoint(idx, mt, dsdt);
			}

			g_dsdt->Write();

			delete im;
		}
		
		IslamModel2015 *im = new IslamModel2015;
		im->Configure(IslamModel2015::mFull);
		im->Init();
		im->Print();
		
		TGraph *g_chi_om_re = new TGraph(); g_chi_om_re->SetName("g_chi_om_re");
		TGraph *g_chi_om_im = new TGraph(); g_chi_om_im->SetName("g_chi_om_im");

		for (double b = 0.; b <= 10.; b += 0.1)
		{
			TComplex chi_om = im->Chi_omega(b);

			int idx = g_chi_om_re->GetN();
			g_chi_om_re->SetPoint(idx, b, chi_om.Re());
			g_chi_om_im->SetPoint(idx, b, chi_om.Im());
		}

		g_chi_om_re->Write();
		g_chi_om_im->Write();

		TGraph *g_F_cal = new TGraph(); g_F_cal->SetName("g_F_cal");

		for (double qt = 0.; qt <= 10.; qt += 0.1)
		{
			TComplex F_cal = im->F_cal(qt);

			int idx = g_F_cal->GetN();
			g_F_cal->SetPoint(idx, qt, F_cal);
		}

		g_F_cal->Write();

		delete f_out;
	}
}

//----------------------------------------------------------------------------------------------------

int main()
{
	TFile *f_out = new TFile("debugger.root", "recreate");

	// initialise constants etc.
	//Constants::Init(53, Constants::mPP);
	Constants::Init(8000, Constants::mPP);
	//Constants::Init(13000, Constants::mPP);
	
	coulomb->precision = 1E-5;
	coulomb->mode = CoulombInterference::mKL;

	coulomb->InitIFunctionInterpolator(20., 300);

	cnts->Print();
	coulomb->Print();

	ModelFactory mf;
	model = mf.MakeInstance("block [06]");
	//model = mf.MakeInstance("bourrely [03]", false);
	//model = mf.MakeInstance("dl [13]");
	//model = mf.MakeInstance("godizov [14]", false);
	//model = mf.MakeInstance("islam (hp) [06,09]");
	//model = mf.MakeInstance("islam (lxg) [06,09]");
	//model = mf.MakeInstance("jenkovszky [11]");
	//model = mf.MakeInstance("petrov (2p) [02]");
	//model = mf.MakeInstance("petrov (3p) [02]");
	model->Print();

	// dsdt
	TGraph *g_dsdt = new TGraph();
	g_dsdt->SetName("g_dsdt");
	for (double mt = 1e-4; mt < 1; mt += 0.01)
	{
		double dsdt = cnts->sig_fac * coulomb->Amp(-mt).Rho2();
		//printf("t = %E: dsdt = %E\n", mt, dsdt);

		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, mt, dsdt);
	}
	g_dsdt->Write();

/*
	// profile
	TGraph *g_prof = new TGraph();
	g_prof->SetName("g_prof");
	for (double b = 0.; b < 8.; b += 0.02)
	{
		double prof = model->Prf(b).Rho();
		printf("b = %E: |prof| = %E\n", b, prof);

		int idx = g_prof->GetN();
		g_prof->SetPoint(idx, b, prof);
	}
	g_prof->Write();
*/

	// test A term
	/*
	for (double mt = 1e-4; mt < 0.2; mt += 0.01)
	{
		double r = coulomb->A_term(-mt);
		printf("%.4f, %.3E\n", mt, r);
	}
	*/

	// test I term
	/*
	for (double mt = 1e-4; mt < 0.21; mt += 0.1)
		for (double mt2 = 1e-4; mt2 < 0.21; mt2 += 0.1)
		{
			if (mt == mt2)
				continue;

			double r = coulomb->I_integral(-mt, -mt2);
			printf("%.4f, %.4f | %.3E\n", mt, mt2, r);
		}
	*/

	// test B term
	/*
	for (double mt = 1e-4; mt < 0.2; mt += 0.01)
	{
		printf("* mt = %E\n", mt);
		TComplex r = coulomb->B_term(-mt);
		printf(" => %.3E, %.3E\n", r.Re(), r.Im());
	}
	*/

	TGraph *g_B_re = new TGraph(); g_B_re->SetName("g_B_re");
	TGraph *g_B_im = new TGraph(); g_B_im->SetName("g_B_im");
	for (double mt = 1e-4; mt < 1; mt += 0.01)
	{
		TComplex B = coulomb->B_term(-mt);

		int idx = g_B_re->GetN();
		g_B_re->SetPoint(idx, mt, B.Re());
		g_B_im->SetPoint(idx, mt, B.Im());
	}
	g_B_re->Write();
	g_B_im->Write();

	// C, R, Z graphs
/*
	TGraph *g_C = new TGraph(); g_C->SetName("g_C"); g_C->SetLineColor(2);
	TGraph *g_R = new TGraph(); g_R->SetName("g_R"); g_R->SetLineColor(4);
	TGraph *g_Z = new TGraph(); g_Z->SetName("g_Z"); g_Z->SetLineColor(6);
	for (double mt = 4e-3; mt < 10.; mt += 0.03)
	{
		printf("* mt = %E\n", mt);
		double C = coulomb->C(-mt);
		double R = coulomb->R(-mt);
		double Z = coulomb->Z(-mt);

		int idx = g_C->GetN();
		g_C->SetPoint(idx, mt, C);
		g_R->SetPoint(idx, mt, R);
		g_Z->SetPoint(idx, mt, Z);
	}

	g_C->Write();
	g_R->Write();
	g_Z->Write();
*/

	/*
	// s-dependent quantities
	model->ForcePresampling(false);

	for (unsigned int mi = 0; mi < 2; mi++)
	{
		Constants::ParticleMode pMode = (mi == 0) ? Constants::mPP : Constants::mAPP;

		gDirectory = f_out->mkdir((mi == 0) ? "pp" : "app");

		TGraph *g_si_tot = new TGraph(); g_si_tot->SetName("g_si_tot");
		TGraph *g_rho = new TGraph(); g_rho->SetName("g_rho");
		TGraph *g_B0 = new TGraph(); g_B0->SetName("g_B0");
		for (double W = 10; W <= 1e4; W *= 1.5)
		{
			Constants::Init(W, pMode);
			model->Init();
	
			double ep = 1E-5;
			TComplex amp0 = model->Amp(0.);
			TComplex amp_ep = model->Amp(-ep);
	
			double si_tot = 4.*cnts->pi*cnts->sq_hbarc/cnts->p_cms/cnts->sqrt_s * amp0.Im();
			double rho = (amp0.Im() != 0.) ? amp0.Re() / amp0.Im() : 0.;
			double B0 = ( log(amp0.Rho2()) - log(amp_ep.Rho2()) ) / ep;
	
			int idx = g_si_tot->GetN();
			g_si_tot->SetPoint(idx, W, si_tot);
			g_rho->SetPoint(idx, W, rho);
			g_B0->SetPoint(idx, W, B0);

			printf("W = %.1f, si_tot = %.2f\n", W, si_tot);
		}
	
		g_si_tot->Write();
		g_rho->Write();
		g_B0->Write();
	}
	*/

	delete f_out;
	return 0;
}
