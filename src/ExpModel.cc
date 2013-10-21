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

#include "interface/ExpModel.h"
#include "interface/Constants.h"

using namespace Elegent;

//----------------------------------------------------------------------------------------------------

ExpModel::ExpModel()
{
	fullLabel.name = "exponential"; shortLabel.name = "exp";
}

//----------------------------------------------------------------------------------------------------

void ExpModel::Configure()
{
	fullLabel.variant = ""; shortLabel.variant = "";
	fullLabel.version = ""; shortLabel.version = "";
	fullLabel.mode = ""; shortLabel.mode = "";
}

//----------------------------------------------------------------------------------------------------

void ExpModel::Init()
{
	a = 2E9;
	b1 = 10.;
	b2 = 0.;
	p0 = M_PI/2.;
	p1 = 0.;
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
	printf(">> ExpModel::Prf > not implemented.\n");
	return 0.;
}

//----------------------------------------------------------------------------------------------------

TComplex ExpModel::Amp(double t) const
{
	return a * TComplex::Exp( b1*t + b2*t*t + i*(p0 + p1*t) );
}
