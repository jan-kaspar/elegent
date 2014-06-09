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

#ifndef _elegent_math_
#define _elegent_math_

#include <TComplex.h>

#include <gsl/gsl_integration.h>

namespace Elegent
{

typedef double (* RealFunction)(double x, double *par, const void *obj);
typedef TComplex (* ComplexFunction)(double x, double *par, const void *obj);

// TODO: describe
double RealIntegrate(RealFunction fcn, double *par, const void *object, double from, double to, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space);

// TODO: describe
TComplex ComplexIntegrate(ComplexFunction fcn, double *par, const void *object, double from, double to, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space);

} // namespace

#endif
