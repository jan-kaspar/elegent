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
#include "TGraph.h"

using namespace Elegent;
using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// initialize calculation engine: for proton-proton scattering at sqrt(s) = 13 TeV
	Constants::Init(2*6500., cnts->mPP);
    cnts->Print();

	// initialize hadronic model: e.g. model by Petrov et al. (with 3 Pomerons)
	PPPModel *ppp3 = new PPPModel();
	ppp3->Configure(PPPModel::v3P);
	ppp3->Init();
	model = ppp3;
	model->Print();

	// initialize CNI calculation
	coulomb->mode = CoulombInterference::mKL;	// formula by Kundrat-Lokajicek
	coulomb->ffType = coulomb->ffPuckett;
	coulomb->precision = 1E-3;
	coulomb->Print();

	// sample differential cross-section
	TFile *f_out = TFile::Open("example1.root", "recreate");

	TGraph *g_dsdt = new TGraph();
	g_dsdt->SetName("g_dsdt");

	for (double mt = 1E-5; mt < 1.; mt *= 1.1)
	{
		const double t = -mt;	// four-momentum transfer squared, in GeV^2

		const TComplex A = coulomb->Amp(t);		// scattering amplitude

		const double dsdt = cnts->sig_fac * A.Rho2();	// differential cross-section

		g_dsdt->SetPoint(g_dsdt->GetN(), mt, dsdt);
	}

	g_dsdt->Write();

	delete f_out;

	return 0;
}
