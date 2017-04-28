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

#ifndef _elegent_i_function_interpolator_
#define _elegent_i_function_interpolator_

#include <gsl/gsl_spline2d.h>

namespace Elegent
{

class CoulombInterference;

/**
 * Class to interpolate \f$I(t, t')\f$ pre-evaluated on grid of t vs. t'.
 **/
class IFunctionInterpolator
{
	private:
		/// maximum of |t| for which interpolation works
		double mt_max;
	
		/// interpolation objects
		double *mt_values;
		double *I_values_on_grid;
		gsl_spline2d *spline;
		gsl_interp_accel *xacc, *yacc;

	public:
		/// constructor evaluates \f$I(t, t')\f$ on a grid
		IFunctionInterpolator(CoulombInterference *coulomb, double _mt_max, unsigned int n_grid_values, double t_diff_min=1E-5, bool bicubic=true);

		/// destructor calls Release if previously initialised
		~IFunctionInterpolator();

		/// interpolate \f$I(t, t')\f$ from pre-evaluated grid
		/// \param t in GeV^2, negative
		/// \param tp in GeV^2, negative
		double Eval(double t, double tp) const;
};

} // namespace

#endif
