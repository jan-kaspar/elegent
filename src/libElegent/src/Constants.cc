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

#include "interface/Config.h"
#include "interface/Constants.h"

// TODO: remove to better place
#include <gsl/gsl_errno.h>

namespace Elegent
{

//----------------------------------------------------------------------------------------------------

TComplex i(0, 1);
Constants *cnts = NULL;

// physics constants
double Constants::alpha =			7.297E-3;
double Constants::proton_mass =		0.938271;
double Constants::neutron_mass =	0.939565;
double Constants::hbarc =			0.197326;
double Constants::sq_hbarc =		0.389379;
double Constants::M = proton_mass;
double Constants::M_sq= proton_mass * proton_mass;
double Constants::kappa =			1.793;

// mathematical constants
double Constants::pi =				M_PI;
double Constants::gamma =			0.577215;

//----------------------------------------------------------------------------------------------------

void Constants::Configure(double W, Constants::ParticleMode mode)
{
	sqrt_s = W;
	s = sqrt_s*sqrt_s;
	ln_s = log(s);
	p_cms = sqrt(s /4. - proton_mass * proton_mass);
	sig_fac = sq_hbarc * pi / (s * p_cms * p_cms);
	t_min = 4.*proton_mass * proton_mass - s;
	
	pMode = mode;
}

//----------------------------------------------------------------------------------------------------

void Constants::Init(double W, Constants::ParticleMode mode)
{
	cnts = new Constants(W, mode);

	// TODO: find a better place
	// 	maybe make a sort of Elegent master class, running all-common inits, 
	//	providing also model factory, lists of available models etc.
	gsl_set_error_handler_off();
}

//----------------------------------------------------------------------------------------------------

void Constants::Print()
{
	printf(">> Constants::Print\n");

	// TODO: find a better place
	// 	maybe make a sort of Elegent master class, running all-common inits, 
	//	providing also model factory, lists of available models etc.
	printf("\tElegent version: " Elegent_VERSION "\n");

	printf("\talpha = %E\n", alpha);
	printf("\tproton_mass = %E\n", proton_mass);
	printf("\tneutron_mass = %E\n", neutron_mass);
	printf("\thbarc = %E\n", hbarc);
	printf("\tsq_hbarc = %E\n", sq_hbarc);
	printf("\tM = %E\n", M);
	printf("\tM_sq = %E\n", M_sq);
	printf("\tpi = %E\n", pi);
	printf("\tgamma = %E\n", gamma);
	printf("\tsqrt_s = %E\n", sqrt_s);
	printf("\ts = %E\n", s);
	printf("\tln_s = %E\n", ln_s);
	printf("\tp_cms = %E\n", p_cms);
	printf("\tsig_fac = %E\n", sig_fac);
	printf("\tt_min = %E\n", t_min);
	printf("\tpMode = %s\n", (pMode == mPP) ? "pp" : "app");
}

} // namespace
