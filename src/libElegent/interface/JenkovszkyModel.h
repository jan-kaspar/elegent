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

#ifndef _elegent_jenkovszky_model_
#define _elegent_jenkovszky_model_

#include "Model.h"
#include "interface/Math.h"

namespace Elegent
{

/**
 * The model of Jenkovszky et al.
 * 
 * References:
 *	[1] L. L. JENKOVSZKY, A. I. LENGYEL, and D. I. LONTKOVSKYI, Int. J. Mod. Phys. A 26 (2011) 4755. DOI: 10.1142/S0217751X11054760 
 **/
class JenkovszkyModel : public Model
{
	public:
		JenkovszkyModel();
		~JenkovszkyModel();

		void Configure();
		virtual void Init();

		virtual void Print() const;

		virtual TComplex Amp(double t) const;
		virtual TComplex Prf(double b) const;

	protected:
		/// pomeron parameters
		double a_P, b_P, de_P, al1_P, ep_P, s_P;

		/// odderon parameters
		double a_O, b_O, de_O, al1_O, s_O;

		/// omega parameters
		double a_om, b_om, s_om, al0_om, al1_om;

		/// f parameters
		double a_f, b_f, s_f, al0_f, al1_f;
		
		/// integration parameters for profile-funcion calculation
		double precision_t, upper_bound_t;

		bool integ_workspace_initialized;
		unsigned long integ_workspace_size;
		gsl_integration_workspace *integ_workspace;
		
		static TComplex Amp_J0(double t, double *par, const void *obj);
};

} // namespace

#endif
