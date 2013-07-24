/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#include "interface/InterpolationModel.h"

#include "TGraph.h"

namespace Elegent
{

//#define DEBUG 1

//----------------------------------------------------------------------------------------------------

InterpolationModel::~InterpolationModel()
{
}

//----------------------------------------------------------------------------------------------------

void InterpolationModel::Print() const
{
	printf(">> InterpolationModel::Print\n");
	printf("\tN = %u, t_min = %.3E, t_max = %.3E\n", N, t_min, t_max);
}

//----------------------------------------------------------------------------------------------------

TComplex InterpolationModel::Prf(double) const
{
	// this function is not planned to be used
	printf(">> InterpolationModel::Prf > not implemented.\n");
	return 0.;
}


} // namespace
