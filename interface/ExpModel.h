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

#ifndef _elegent_exp_model_
#define _elegent_exp_model_

#include "Model.h"

namespace Elegent
{

/**
 * Exponential model of p-p and p-anti p elastic scattering.
 *
 * amplitude(t) = a * exp( b1*t + b2*t*t + i*(p0 + p1*t) )
 * For reference pourposes only.
 **/
class ExpModel : public Model
{
	public:
		ExpModel();
		
		void Configure();
		virtual void Init();
		virtual void Print() const;
		virtual TComplex Amp(double t) const;
		virtual TComplex Prf(double b) const;

	public:
		double a, b1, b2, p0, p1;
};

} // namespace

#endif
