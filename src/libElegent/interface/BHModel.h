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

#ifndef _elegent_bh_model_
#define _elegent_bh_model_

#include "Model.h"
#include "Math.h"

namespace Elegent
{

/**
 * Block-Halzen model of p-p and p-anti p elastic scattering.
 *
 * References:
 *	[1] BLOCK, M. M., GREGORES, E. M., HALZEN, F. and PANCHERI, G., Phys. Rev. D60 (1999) 054024
 *	[2] BLOCK, M. M., Phys. Rept. 436 (2006) 71-215
 **/
class BHModel : public Model
{
	public:
		BHModel();
		~BHModel();
		
		void Configure();
		virtual void Init();
		virtual void Print() const;
		virtual TComplex Amp(double t) const;
		virtual TComplex Prf(double b) const;

	protected:
		/// common parameters
		double m0, s0, al_s, Sigma_gg; 

		/// parameters for sigma_gg
		double Cp_gg, epsilon, Ng; 
		std::vector<double> a, b;
		
		/// parameters for sigma_gg
		double C_qg_log; 
		
		/// parameters for sigma_qg
		double C, C_even_regge; 
		
		/// parameters for sigma_qq
		double C_odd; 

		/// effective areas
		double mu_gg, mu_qg, mu_qq, mu_odd;

		/// the integral cross-sections
		TComplex sigma_gg, sigma_qq, sigma_qg, sigma_odd;
		
		/// integration parameters
		double precision, upper_bound;

		bool integ_workspace_initialized;
		unsigned long integ_workspace_size;
		gsl_integration_workspace *integ_workspace;

		/// profile defined by Eq. (B2) in [1]
		double W(double b, double mi) const;

		/// sum in Eq. (B5) in [1]
		TComplex Sum(double s) const;

		/// the full eikonal: Eq. (1) without the leading factor i and Eq. (12)
		TComplex chi_without_i(double b) const;

		TComplex prf0(double b) const;
		static TComplex prf0_J0(double b, double *par, const void *vobj);
};

} // namespace

#endif
