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

#ifndef _elegent_godizov_model_
#define _elegent_godizov_model_

#include "Model.h"
#include "Math.h"

namespace Elegent
{

/**
 * \brief TODO
 * References:
 *	[1] arXiv:1404.2851v2
 **/
class GodizovModel : public Model
{
	public:
		GodizovModel();
		~GodizovModel();
		
		void Configure(bool _presampled = true);

		virtual void Init();

		virtual void Print() const;

		virtual TComplex Amp(double t) const;

		/// b in fm
		virtual TComplex Prf(double b) const;

	protected:
		double De;		///< al_P(0) - 1
		double ta_a;
		double Ga_P0;
		double ta_g;
		double s0;

		/// flag whether the profile function is presampled
		bool presampled;

		/// Eq. (2) in [1]
		TComplex delta_t(double t) const;

		static TComplex delta_t_J0(double t, double *par, const void *obj);
		
		/// bottom relation from Eq. (1) in [1]
		TComplex delta_b(double b) const;

		/// profile function, b in GeV^-1, see Eq. (1) in [1]
		TComplex prf0(double b) const;

		static TComplex prf_J0(double b, double *par, const void *obj);

		/// integration variables
		double upper_bound_t, precision_t;
		double upper_bound_b, precision_b;

		bool integ_workspace_initialized;
		unsigned long integ_workspace_size_b;
		gsl_integration_workspace *integ_workspace_b;
		unsigned long integ_workspace_size_t;
		gsl_integration_workspace *integ_workspace_t;

		/// the sampling-step size
		double prf0_sample_db;

		/// the number of sampled points
		unsigned int prf0_sample_N;

		/// the sampled real and imaginary values of prf0(b) 
		std::vector<double> prf0_sample_re, prf0_sample_im;

		/// samples the prf0 function
		void Prf0SampleBuild(unsigned int samples);

		/// interpolates (linearly) the sample at point b
		TComplex Prf0SampleEval(double b) const;
};

} // namespace

#endif
