// HEADER

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
	printf(">> InterpolationModel::Prf > Not implemented.\n");
	return 0;
}


} // namespace
