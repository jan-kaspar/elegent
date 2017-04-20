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

#ifndef _elegent_interpolation_model_
#define _elegent_interpolation_model_

#include "Model.h"

#include <vector>

namespace Elegent
{

/**
 * Model that interpolates stored amplitude points.
 **/
class InterpolationModel : public Model
{
	public:
		/// number of points (samples)
		unsigned int N;

		/// the lower boundary
		double t_min;

		/// the upper boundary
		double t_max;

		/// the t interval between two adjacent (equidistant) points
		double dt;

		/// amplitude samples
		std::vector<TComplex> amp_data;

	public:
		InterpolationModel(unsigned int _N, double _t_min, double _t_max);

		virtual ~InterpolationModel();

		void Configure();
		
		virtual void Init() {}

		virtual void Print() const;

		/// amplitude, t in GeV^-2, t < 0
		virtual TComplex Amp(double t) const
		{
			double f = (t - t_min) / dt;
			unsigned int idx = (unsigned int) f;
			
			if (idx + 1 > N - 1)
			{
				if (fabs(t - t_max) < 1E-10)
					return amp_data[N-1];
				return TComplex(0, 0);
			}
		
			f -= idx;
		
			return amp_data[idx] + (amp_data[idx+1] - amp_data[idx]) * f;
		}

		virtual TComplex Prf(double b) const;

		/// returns the value of t corresponding to the point with index `idx' (between 0 and N-1 inclusively)
		inline double GetT(unsigned int idx) const
		{
			return t_min + dt * idx;
		}

		inline void SetPoint(unsigned int idx, const TComplex &v)
		{
			amp_data[idx] = v;
		}
};

} // namespace

#endif
