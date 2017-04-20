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

#ifndef _elegent_islam_model_2015_
#define _elegent_islam_model_2015_

#include "Model.h"
#include "Math.h"

namespace Elegent
{


/**
 * Islam model of p-p and p-anti p elastic scattering.
 * 
 * References:
 *	[1] EDS'15 paper draft
 *	[2] private communication
 **/
class IslamModel2015 : public Model
{
	public:
		/// mode of the model
		enum ModeType
		{
			mDiff,			///< diffraction amplitude
			mOmega,			///< omega-exchange amplitude
			mQuark,			///< quark-quark amplitude
			mPolarization,	///< polarisation amplitude
			mFull			///< diffraction, core and quark-quark amplitude
		} mode;

		IslamModel2015();
		~IslamModel2015();

		void Configure(ModeType _m);

		virtual void Init();

		virtual void Print() const;

		virtual TComplex Amp(double t) const;
		virtual TComplex Prf(double b) const;
		
	public:
		/// diffraction amplitude
		double R0, R1, a0, a1;
		TComplex R, a;
		double eta0, c0, sigma;
		TComplex g_s;

		TComplex GammaDPlus(double b) const;
		static TComplex GammaDPlus_J0(double b, double *par, const void *vobj);
		TComplex T_diff(double t) const;

		/// screening due to diffraction
		double lambda0, d0, s0;
		TComplex fac_scr_diff;	///< screening factor due to diffraction
		
		/// omega-exchange amplitude
		double gamma_hat, theta_hat, beta, m_omega;

		TComplex Chi_omega(double b) const;
		static TComplex GammaOmega_J0(double b, double *par, const void *vobj);
		TComplex T_omega(double t) const;

		/// screening due to omega exchange
		double b_tilde;
		
		TComplex fac_scr_omega;	///< screening factor due to omega exange
	
		/// quark confinement parameters
		double m0, M;
	
		/// "low-x gluons" amplitude
		double gamma_gg, lambda, m, mu;

		double I_integral(double qt, double al) const;
		static double F_cal_integ(double x, double *par, const void *vobj);
		double F_cal(double qt) const;

		TComplex T_quark(double t) const;

		/// polarisation amplitude
		double ht, wd, sf;

		TComplex T_polarization(double t) const;
	
		/// integration variables
		double precision_b, precision_t, upper_bound_b, upper_bound_t;

		bool integ_workspace_initialized;
		unsigned long integ_workspace_size_b;
		gsl_integration_workspace *integ_workspace_b;
		unsigned long integ_workspace_size_t;
		gsl_integration_workspace *integ_workspace_t;

		/// profile funcion methods
		static TComplex Amp_J0(double t, double *par, const void *vobj);
};

} // namespace

#endif
