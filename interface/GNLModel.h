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

#ifndef _elegent_gnl_model_
#define _elegent_gnl_model_

#include "Model.h"
#include "Math.h"

namespace Elegent
{

/**
 * Model of elastic pp scattering by Gauron, Nicolescu, Leader et al.
 * 
 * References:
 *	[1] arXiv:hep-ph/0607089v4
 **/
class GNLModel : public Model
{
	public:
		GNLModel();
		~GNLModel();
		
		void Configure();

		virtual void Init() override;

		virtual void Print() const override;

		virtual TComplex Amp(double t) const override;

		/// b in fm
		virtual TComplex Prf(double b) const override;

	protected:
		/// parameters of F_H_pl
		double H_1, b_1_pl, H_2, b_2_pl, H_3, b_3_pl, K_pl;

		/// parameters of Odderon
		double O_1, b_1_mi, O_2, b_2_mi, O_3, b_3_mi, K_mi;

		struct Singularity
		{
			double al_0, C, be, al_p;
		};

		Singularity s_P, s_PP, s_O, s_OP, s_R_pl, s_R_mi, s_RP_pl, s_RP_mi;

		/// integration variables
		double upper_bound_t, precision_t;

		bool integ_workspace_initialized;
		unsigned long integ_workspace_size_t;
		gsl_integration_workspace *integ_workspace_t;
};

} // namespace

#endif
