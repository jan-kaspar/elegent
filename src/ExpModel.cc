// HEADER

#include "interface/ExpModel.h"
#include "interface/Constants.h"

using namespace Elegent;

//----------------------------------------------------------------------------------------------------

ExpModel::ExpModel() : Model("exp", "Exp", "exponential (reference)"),
  a(2E9), b1(10.), b2(0), p0(M_PI/2), p1(0)
{
}

//----------------------------------------------------------------------------------------------------

void ExpModel::Print() const
{
  printf(">> ExpModel::Print\n");
  printf("\ta=%E\n", a);
  printf("\tb1=%E, b2=%E\n", b1, b2);
  printf("\tp0=%E, p1=%E\n", p0, p1);
}

//----------------------------------------------------------------------------------------------------

TComplex ExpModel::Prf(double ) const
{
  // this function is not planned to be used
  printf(">> ExpModel::Prf > Not implemented\n");
  return 0;
}

//----------------------------------------------------------------------------------------------------

TComplex ExpModel::Amp(double t) const
{
  return  a * TComplex::Exp( b1*t + b2*t*t + i*(p0 + p1*t) );
}

// vim: ft=cmscpp
