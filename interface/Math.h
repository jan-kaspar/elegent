/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_math_
#define _elegent_math_

#include <TComplex.h>

namespace Elegent
{

// ================================== INTEGRATION ROUTINES =================================================

// TODO: description
double DoubleInt(const void *obj, double (*fcn)(double*, double*, const void*), double a, double b, double *params = NULL, double epsilon = 1E-9);
TComplex CmplxInt(const void *obj, TComplex (*fcn)(double*, double*, const void *), double a, double b, double *params = NULL, double epsilon = 1E-9);

} // namespace

#endif
