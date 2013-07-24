/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_math_
#define _elegent_math_

#include <TComplex.h>

namespace Elegent
{

// =================================== DERIVATION ROUTINES =================================================

// TODO: needed?
TComplex CmplxDeriv(TComplex (*fcn)(Double_t*, Double_t*), Double_t x, Double_t *params = NULL, Double_t h = 1E-5);

// ================================== INTEGRATION ROUTINES =================================================

Double_t DoubleInt(const void *obj, Double_t (*fcn)(Double_t*, Double_t*, const void*), Double_t a, Double_t b, Double_t *params = NULL, Double_t epsilon = 1E-9);
TComplex CmplxInt(const void *obj, TComplex (*fcn)(Double_t*, Double_t*, const void *), Double_t a, Double_t b, Double_t *params = NULL, Double_t epsilon = 1E-9);

} // namespace

#endif
