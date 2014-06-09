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


namespace Elegent
{

//----------------------------------------------------------------------------------------------------

// TODO: better names
// TODO: describe

struct RealIntegPar
{
	RealFunction fcn;
	double *parameters;
	const void *object;

	RealIntegPar(RealFunction _f, double *_p, const void *_o) : fcn(_f), parameters(_p), object(_o) {}
};

//----------------------------------------------------------------------------------------------------

double RealIntegFcn(double x, void *vpar)
{
	RealIntegPar *par = (RealIntegPar *) vpar;
	return par->fcn(x, par->parameters, par->object);
}

//----------------------------------------------------------------------------------------------------

double RealIntegrate(RealFunction fcn, double *par, const void *object, double from, double to, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space)
{
	// prepare structures
	RealIntegPar ocip(fcn, par, object);

	gsl_function F;
  	F.function = RealIntegFcn;
  	F.params = &ocip;

	// real part
	double result, error;
	gsl_integration_qag(&F, from, to, 0., rel_err, work_space_size, GSL_INTEG_GAUSS61, work_space, &result, &error);

	return result;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

struct OneCompIntegPar
{
	ComplexFunction fcn;
	double *parameters;
	const void *object;
	enum Part { pUndefined, pReal, pImaginary } part;

	OneCompIntegPar(ComplexFunction _f, double *_p, const void *_o, Part _pa = pUndefined) : fcn(_f), parameters(_p), object(_o), part(_pa) {}
};

//----------------------------------------------------------------------------------------------------

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

TComplex ComplexIntegrate(ComplexFunction fcn, double *par, const void *object, double from, double to, double rel_err,
	unsigned long work_space_size, gsl_integration_workspace *work_space)
{
	// prepare structures
	OneCompIntegPar ocip(fcn, par, object);

	gsl_function F;
  	F.function = OneCompIntegFcn;
  	F.params = &ocip;

	// real part
	ocip.part = OneCompIntegPar::pReal;
	double result_re, error_re;
	gsl_integration_qag(&F, from, to, 0., rel_err, work_space_size, GSL_INTEG_GAUSS61, work_space, &result_re, &error_re);

	// imaginary part
	ocip.part = OneCompIntegPar::pImaginary;
	double result_im, error_im;
	gsl_integration_qag(&F, from, to, 0., rel_err, work_space_size, GSL_INTEG_GAUSS61, work_space, &result_im, &error_im);

	// result
	return TComplex(result_re, result_im);
}
	

} // namespace
