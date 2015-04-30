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

#include "interface/Math.h"

#include <gsl/gsl_errno.h>


namespace Elegent
{

/// false by default
bool abortOnIntegrationError = false;

//----------------------------------------------------------------------------------------------------

/// Set of parameters for integration of a real function.
/// \ingroup Math
struct RealIntegPar
{
	RealFunction fcn;		///< function to integrate
	double *parameters;		///< additional numerical parameters
	const void *object;		///< additional pointer to an object

	RealIntegPar(RealFunction _f, double *_p, const void *_o) : fcn(_f), parameters(_p), object(_o) {}
};

//----------------------------------------------------------------------------------------------------

/// Evaluates a given real function with a given parameter set.
/// \ingroup Math
double RealIntegFcn(double x, void *vpar)
{
	RealIntegPar *par = (RealIntegPar *) vpar;
	return par->fcn(x, par->parameters, par->object);
}

//----------------------------------------------------------------------------------------------------

double RealIntegrate(RealFunction fcn, double *par, const void *object,
	double from, double to,
	double abs_err, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space, const char *errorLabel)
{
	// prepare structures
	RealIntegPar ocip(fcn, par, object);

	gsl_function F;
  	F.function = RealIntegFcn;
  	F.params = &ocip;

	// real part
	double result, unc;
	int status = gsl_integration_qag(&F, from, to, abs_err, rel_err, work_space_size, GSL_INTEG_GAUSS61, work_space, &result, &unc);
	if (status != 0) 
	{
		printf("WARNING in %s > Integration failed: %s.\n", errorLabel, gsl_strerror(status));
		if (abortOnIntegrationError)
			abort();
	}

	return result;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

/// Set of parameters for integration of a part of a complex function.
/// \ingroup Math
struct OneCompIntegPar
{
	ComplexFunction fcn;	///< function to integrate
	double *parameters; 	///< additional numerical parameters
	const void *object; 	///< additional pointer to an object

	/// part to integrate
	enum Part { pUndefined, pReal, pImaginary } part;

	OneCompIntegPar(ComplexFunction _f, double *_p, const void *_o, Part _pa = pUndefined) : fcn(_f), parameters(_p), object(_o), part(_pa) {}
};

//----------------------------------------------------------------------------------------------------

/// Evaluates a given part of a give complex function with a given parameter set.
/// \ingroup Math
double OneCompIntegFcn(double x, void *par)
{
	OneCompIntegPar *ocip = (OneCompIntegPar *) par;
	TComplex c = ocip->fcn(x, ocip->parameters, ocip->object);

	if (ocip->part == OneCompIntegPar::pReal)
		return c.Re();
	if (ocip->part == OneCompIntegPar::pImaginary)
		return c.Im();

	return 0.;
}

//----------------------------------------------------------------------------------------------------

TComplex ComplexIntegrate(ComplexFunction fcn, double *par, const void *object,
	double from, double to,
	double abs_err, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space, const char *errorLabel)
{
	// prepare structures
	OneCompIntegPar ocip(fcn, par, object);

	gsl_function F;
  	F.function = OneCompIntegFcn;
  	F.params = &ocip;

	// real part
	ocip.part = OneCompIntegPar::pReal;
	double result_re, unc_re;
	int status_re = gsl_integration_qag(&F, from, to, abs_err, rel_err, work_space_size, GSL_INTEG_GAUSS61, work_space, &result_re, &unc_re);
	if (status_re != 0) 
	{
		printf("WARNING in %s > Real integration failed: %s.\n", errorLabel, gsl_strerror(status_re));
		if (abortOnIntegrationError)
			abort();
	}

	// imaginary part
	ocip.part = OneCompIntegPar::pImaginary;
	double result_im, unc_im;
	int status_im = gsl_integration_qag(&F, from, to, abs_err, rel_err, work_space_size, GSL_INTEG_GAUSS61, work_space, &result_im, &unc_im);
	if (status_im != 0) 
	{
		printf("WARNING in %s > Imaginary integration failed: %s.\n", errorLabel, gsl_strerror(status_im));
		if (abortOnIntegrationError)
			abort();
	}

	// result
	return TComplex(result_re, result_im);
}
	

} // namespace
