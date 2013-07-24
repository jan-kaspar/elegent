// HEADER

#ifndef _elegent_interpolation_model_
#define _elegent_interpolation_model_

#include "Model.h"

#include <vector>

class TGraph;

namespace Elegent
{

/**
 * \ingroup Elegent
 * \brief Model that interpolates stored amplitude points.
 **/
class InterpolationModel : public Model
{
	public:
		/// TODO
		unsigned int N;

		/// the lower boundary
		double t_min;

		/// the upper boundary
		double t_max;

		/// TODO
		double dt;

		/// TODO
		std::vector<TComplex> amp_data;

	public:
		InterpolationModel(unsigned int _N, double _t_min, double _t_max) :
			N(_N), t_min(_t_min), t_max(_t_max), dt( (t_max - t_min) / (N-1) ), amp_data(N)
		{
		}

		virtual ~InterpolationModel();

		virtual std::string GetModeString() const
			{ return "basic"; }

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

		/// TODO
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
