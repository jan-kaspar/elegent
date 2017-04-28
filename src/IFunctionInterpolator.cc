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

#include "interface/IFunctionInterpolator.h"

#include "interface/CoulombInterference.h"

namespace Elegent
{

IFunctionInterpolator::IFunctionInterpolator(CoulombInterference *coulomb, double _mt_max, unsigned int n_grid_values, double t_diff_min, bool bicubic)
{
	mt_max = _mt_max;

	// define grid
	mt_values = new double[n_grid_values];
	for (unsigned int i = 0; i < n_grid_values; i++)
	{
		const double p = 2.;
		mt_values[i] = mt_max * pow(double(i) / (n_grid_values-1), p);
	}

	// create spline object
	const gsl_interp2d_type *interpolation_type = (bicubic) ? gsl_interp2d_bicubic : gsl_interp2d_bilinear;
	spline = gsl_spline2d_alloc(interpolation_type, n_grid_values, n_grid_values);

	// evaluate I_integral on the grid points
	I_values_on_grid = new double[n_grid_values * n_grid_values];
	for (unsigned int i = 0; i < n_grid_values; i++)
	{
		for (unsigned int j = 0; j < n_grid_values; j++)
		{
			double t = -mt_values[i];
			double tp = -mt_values[j];

			double ir = 1.;
			if (fabs(t - tp) > t_diff_min)
			{
				const double I = coulomb->I_function_integration(t, tp);
				const double I0 = -2.*M_PI / fabs(t - tp);
				ir = I/I0;
			}

			double v = ir;
			if (tp < t)
				v = 2. - ir;

			gsl_spline2d_set(spline, I_values_on_grid, i, j, v);
		}
	}

	// initialise interpolation
	gsl_spline2d_init(spline, mt_values, mt_values, I_values_on_grid, n_grid_values, n_grid_values);

	// initialise accelerator objects for look-ups in during the interpolation
	xacc = gsl_interp_accel_alloc();
	yacc = gsl_interp_accel_alloc();
}

//----------------------------------------------------------------------------------------------------

IFunctionInterpolator::~IFunctionInterpolator()
{
	delete[] mt_values;
	delete[] I_values_on_grid;

	gsl_spline2d_free(spline);
	gsl_interp_accel_free(xacc);
	gsl_interp_accel_free(yacc);
}

//----------------------------------------------------------------------------------------------------

double IFunctionInterpolator::Eval(double t, double tp) const
{
	// range check
	if (t > 0. || -t > mt_max || tp > 0. || -tp > mt_max)
		return 0.;

	// do interpolation
	const double v = gsl_spline2d_eval(spline, -t, -tp, xacc, yacc);
	const double ir = (tp < t) ? 2. - v : v;
	const double I0 = -2.*M_PI / fabs(t - tp);

	return ir * I0;
}

} // namespace
