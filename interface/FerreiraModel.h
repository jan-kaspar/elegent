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

#ifndef _elegent_ferreira_model_
#define _elegent_ferreira_model_

#include "Model.h"
#include "Math.h"

namespace Elegent
{

/**
 * Model of elastic pp scattering by Ferreira et al.
 *
 * References:
 *	[1] arXiv:1408.1599v1
 *	[2] A. K. Kohara, E. Ferreira, T. Kadama: Eur. Phys. J. C74 (2014) 11, 3175
 **/
class FerreiraModel : public Model
{
	public:
		FerreiraModel();
		~FerreiraModel();
		
		void Configure();
		virtual void Init();
		virtual void Print() const;
		virtual TComplex Amp(double t) const;
		virtual TComplex Prf(double b_fm) const;

	protected:
		double a0;

		double al_i, al_r;
		double be_i, be_r;
		double ga_i, ga_r;
		double la_i, la_r;

		/// normalisation factor sqrt(s p^2 / pi)
		double N;

		/// Eq. (19) in [1]
		double Psi(double ga, double t) const;

		/// Eq. (14) in [2]
		double R_ggg(double t) const;

		static TComplex Amp_J0(double t, double *par, const void *vobj);

		/// integration variables
		double upper_bound_t, precision_t;

		bool integ_workspace_initialized;
		unsigned long integ_workspace_size_t;
		gsl_integration_workspace *integ_workspace_t;
};

} // namespace

#endif
