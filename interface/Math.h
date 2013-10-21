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

namespace Elegent
{

// ================================== INTEGRATION ROUTINES =================================================

/**
 *\brief Integrates real-valued function `fcn' between points `a' and `b'.
 * `obj' is passed as the third parameter to `fcn'
 * `params' is passed as the second parameter to `fcn'
 * the integration variable is the first element of the array in the first parameter of `fcn'
 * `epsilon' controls the precision of the integration
 *
 * The function has been adapted from ROOT's TF1::Integral.
 **/
double DoubleInt(const void *obj, double (*fcn)(double*, double*, const void*), double a, double b, double *params = NULL, double epsilon = 1E-9);

/**
 *\brief Integrates complex-valued function `fcn' between points `a' and `b'.
 * `obj' is passed as the third parameter to `fcn'
 * `params' is passed as the second parameter to `fcn'
 * the integration variable is the first element of the array in the first parameter of `fcn'
 * `epsilon' controls the precision of the integration
 *
 * The function has been adapted from ROOT's TF1::Integral.
 **/
TComplex CmplxInt(const void *obj, TComplex (*fcn)(double*, double*, const void *), double a, double b, double *params = NULL, double epsilon = 1E-9);

} // namespace

#endif
