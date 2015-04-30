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

/// Abort when integration error occurs?
extern bool abortOnIntegrationError;

/// \brief represents a real function of real variable (\param x)
/// \param par is intended for additional numerical parameters
/// \param obj is inteded for an object pointer
typedef double (* RealFunction)(double x, double *par, const void *obj);

/// \brief integration of a real function of a real variable
/// \param fcn: function to integrate
/// \param par, \param object: additional parameters passed to the integrated function
/// \param from, \param to: integration range
/// \param abs_err, \param rel_err: requested absolute and relative error on the result
/// \param work_space_size, \param work_space: size and reference to the GSL integration workspace
/// \param errorLabel: caption of (possible) error messages
double RealIntegrate(RealFunction fcn, double *par, const void *object,
	double from, double to,
	double abs_err, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space, const char* errorLabel="");

/// \brief represents a complex function of real variable (\param x)
typedef TComplex (* ComplexFunction)(double x, double *par, const void *obj);

/// \brief integration of a complex function of a real variable
TComplex ComplexIntegrate(ComplexFunction fcn, double *par, const void *object,
	double from, double to,
	double abs_err, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space, const char* errorLabel="");

} // namespace

#endif
