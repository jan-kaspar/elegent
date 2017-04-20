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

#include "interface/InterpolationModel.h"

namespace Elegent
{

//#define DEBUG 1

//----------------------------------------------------------------------------------------------------

InterpolationModel::InterpolationModel(unsigned int _N, double _t_min, double _t_max) :
	N(_N), t_min(_t_min), t_max(_t_max), dt( (t_max - t_min) / (N-1) ), amp_data(N)
{
	fullLabel.name = "Interpolation"; shortLabel.name = "itpl";
}

//----------------------------------------------------------------------------------------------------

void InterpolationModel::Configure()
{
}

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
