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

int main()
{
	TFile *f_out = new TFile("debugger.root", "recreate");

	// initialise constants etc.
	Constants::Init(8000, Constants::mPP);
	
	coulomb->precision = 1E-2;

	cnts->Print();
	coulomb->Print();


	ModelFactory mf;
	//model = mf.MakeInstance("petrov (3p) [02]");
	model = mf.MakeInstance("godizov [14]");
	model->Print();

	// dsdt
	TGraph *g_dsdt = new TGraph();
	g_dsdt->SetName("g_dsdt");
	for (double mt = 1e-4; mt < 1.; mt += 0.01)
	{
		double dsdt = cnts->sig_fac * model->Amp(-mt).Rho2();
		printf("t = %E: dsdt = %E\n", mt, dsdt);

		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, mt, dsdt);
	}
	g_dsdt->Write();

	/*
	// profile
	TGraph *g_prof = new TGraph();
	g_prof->SetName("g_prof");
	for (double b = 0.; b < 8.; b += 0.1)
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

	delete f_out;
	return 0;
}
