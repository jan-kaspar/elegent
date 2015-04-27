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

#ifndef _elegent_dl_model_
#define _elegent_dl_model_

#include "Model.h"
#include "Math.h"

namespace Elegent
{

/**
 * \brief Model of elastic pp scattering by Donnachie and Landshoff.
 * References:
 *	[1] Physics Letters B 727 (2013) 500-505
 **/
class DLModel : public Model
{
	public:
		DLModel();
		~DLModel();
		
		void Configure();

		virtual void Init();

		virtual void Print() const;

		virtual TComplex Amp(double t) const;

		/// b in fm
		virtual TComplex Prf(double b) const;

	protected:
		double ep_P, ep_pl, ep_mi;
		double X_P, X_pl, X_mi;
		double al_P_p, al_pl_p, al_mi_p;
		double A, a, b;
		double C, t0;
		double nu;

		/// integration variables
		double upper_bound_t, precision_t;

		bool integ_workspace_initialized;
		unsigned long integ_workspace_size_t;
		gsl_integration_workspace *integ_workspace_t;

		/// Eq. (1c) in [1]
		double F(double t) const;

		/// Eq. (1b) in [1]
		TComplex A_single(double t) const;

		/// Eq. (2b) in [1]
		TComplex A_PP(double t) const;

		/// Eqs. (3a) and (3b) in [1]
		TComplex A_ggg(double t) const;
};

} // namespace

#endif
