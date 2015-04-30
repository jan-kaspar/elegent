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

#ifndef _elegent_constants_
#define _elegent_constants_

#include "TComplex.h"

namespace Elegent
{

/**
 * Set of constants used in Elegent calculations.
 **/
struct Constants
{
	// physics constants
	static double alpha;					///< fine structure constant
	static double proton_mass;				///< GeV
	static double neutron_mass;				///< GeV
	static double hbarc;					///< GeV * fm
	static double sq_hbarc;					///< GeV^2 * mbarn
											///< sigma/mbarn = sigma/GeV^-2 * sq_hbarc
	static double M;						///< abbreviation for proton mass in GeV
	static double M_sq;						///< proton mass squared, GeV^2
	static double kappa;					///< the anomalous magnetic moment of proton
	
	// mathematics constants
	static double pi;						///< pi
	static double gamma;					///< Euler's constant

	// physics data
	double sqrt_s;							///< sqrt_s / GeV
	double s;								///< s / GeV^2
	double ln_s;							///< ln(s / GeV^2)
	double p_cms;							///< particle CMS impuls
	double sig_fac;							///< d sig/dt = sig_fac * |A|^2
	double t_min;							///< negative
	
	/// particle mode
	enum ParticleMode {mPP, mAPP} pMode;
 
	Constants(double W = 0., ParticleMode mode = mPP)
	{
		Configure(W, mode);
	}
	
	/// Configure the constants.
	/// \param W = sqrt(s)
	void Configure(double W, ParticleMode mode);

	/// initilize new instance of Constants and save its pointer to the `cnts' global variable.
	/// \param W = sqrt(s)
	static void Init(double W, ParticleMode mode);

	/// print the actual values
	void Print();
};


extern TComplex i;
extern Constants *cnts;

} // namespace

#endif
