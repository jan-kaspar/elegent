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

/**
 * \defgroup Math
 * Group of utility functions for numerical integration.
 **/

/// Abort when integration error occurs?
/// \ingroup Math
extern bool abortOnIntegrationError;

/// Represents a real function of real variable.
/// \ingroup Math
/// \param x is the integration variable
/// \param par is intended for additional numerical parameters
/// \param obj is inteded for an object pointer
typedef double (* RealFunction)(double x, double *par, const void *obj);

/// Integration of a real function of a real variable.
/// \ingroup Math
double RealIntegrate(
		RealFunction fcn,							///< function to integrate
		double *par,								///< optional numerical parameters
		const void *object,							///< optional object parameter
		double from,								///< lower integration bound
		double to,									///< upper integration bound
		double abs_err,								///< requested absolute error on the result
		double rel_err,								///< requested relative error on the result
		unsigned long work_space_size,				///< size to GSL integration workspace
		gsl_integration_workspace *work_space,		///< pointer to GSL integration workspace
		const char* errorLabel=""					///< caption of (possible) error messages
	);

/// Represents a complex function of real variable.
/// \ingroup Math
typedef TComplex (* ComplexFunction)(double x, double *par, const void *obj);

/// Integration of a complex function of a real variable.
/// \ingroup Math
TComplex ComplexIntegrate(ComplexFunction fcn, double *par, const void *object,
	double from, double to,
	double abs_err, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space, const char* errorLabel="");

} // namespace

#endif
