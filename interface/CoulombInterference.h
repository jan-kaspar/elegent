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

#ifndef _elegent_coulomb_
#define _elegent_coulomb_

#include "Model.h"
#include "Math.h"


namespace Elegent
{

class IFunctionInterpolator;

/**
 * Coulomb hadron interference for elastic scattering.
 *
 * References:
 *  [1] WEST, G. B. and YENNIE, D. R., Phys. Rev. 172 (1968) 1413-1422
 *  [2] KUNDRÁT, V. and LOKAJÍČEK, M., Z. Phys. C63 (1994) 619-630
 *  [3] CAHN, R., Z. Phys. C15 (1982) 23-260
 * TODO: note about wrong sign in WY publication
 * TODO: note about wrong sign in many KL publications
 **/
class CoulombInterference
{
	public:
		/// the mode of coulomb interference
		enum CIMode
		{
			mPC,				///< pure electromagnetic amplitude (Born/OPE approximation) [default]
			mPH,				///< pure hadronic amplitude
			mWY,				///< WY formula
			mSWY,			 	///< simplified WY formula
			mKL,			 	///< (corrected) KL formula
			mCahn				///< Eq. (30) from [3]
		} mode;

		std::string GetModeString() const;

		/// form factor type
		enum FFType
		{
			ffNone,				///< form factor = 1
			ffDipole,			///< dipole form factor, in G_eff, G_E and G_M (G_M(0) = 1)
			ffHofstadter,	 	///< Hofstader et al.: Rev. Mod. Phys. 30 (1958)
			ffBorkowski,		///< Borkowski et al.: Nucl. Phys. B93 (1975)
			ffKelly,			///< Kelly: Phys. Rev. C70 (2004)
			ffArrington,		///< Arrington et al.: Phys. Rev C76 (2007)
			ffPuckett,			///< Puckett et al.: arXiv 1008.0855v1 [default]
			ffPuckettEl		 	///< only electric form-factor of Puckett et al.
		} ffType;

		CoulombInterference();

		~CoulombInterference();

		double GetT()
			{ return T; }

		/// the size of the region around t=t' which is cut off from integration, see B_term method
		double tau;
		
		/// the upper bound of the integration in A_term and B_term is |t|+T
		double T;

		/// precision of the integration
		double precision;

		/// flag whether to use IFunctionInterpolator instead of numerical integration (I_function_integration)
		bool useIFunctionInterpolator = false;

		/// initialises the IFunctionInterpolator
		void InitIFunctionInterpolator(double mt_max, unsigned int n_grid_values);

		/// print the parameters
		void Print() const;

	protected:
		unsigned long integ_workspace_size;
		gsl_integration_workspace *integ_workspace, *integ_workspace2;

		/// optional object for interpolation of \f$I(t, t')\f$
		IFunctionInterpolator *iFunctionInterpolator = NULL;

		/// the integrand of the A term
		static double A_integrand(double tt, double *par, const void *vobj);

		/// the integrand of I(t, t') integral
		static double I_integrand(double phi, double *par, const void *vobj);

		/// the integrand of the B term
		static TComplex B_integrand(double tp, double *par, const void *vobj);

		/// Cahn's B term: expression for integration over phi.
		static TComplex B_cahn_integrand_phi(double phi, double *par, const void *vobj);
		
		/// Cahn's B term: expression for integration over t.
		static TComplex B_cahn_integrand_t(double tp, double *par, const void *vobj);
		
	public:
		/// A: \f$\int_{t_{min}}^0 \log(t'/t) * d/dt(FF^2(t'))\f$.
		/// \param t in GeV^2, negative
		double A_term(double t) const;

		/// evaluates \f$I(t, t')\f$ directly by integration
		double I_function_integration(double t, double tp) const;

		/// evaluates \f$I(t, t')\f$, dependening on useIFunctionInterpolator either by integration or by interpolaion
		double I_function(double t, double tp) const;
		
		/// B: \f${1 / 2\pi} \int_{t_{min}}^0 [ F^N(t') / F^N(t) - 1] I(t, t')\f$.
		/// \param t in GeV^2, negative
		TComplex B_term(double t) const;
		
		/// Cahn's expression for the B term.
		TComplex B_term_cahn(double t) const;

		/// C: the correction for non-vanishing form factors at t_min.
		/// \f$FF^2(t_{min} \log(t/t_{min}))\f$
		/// \param t in GeV^2, negative
		double C_term(double t) const;

		//-------------------- form factors --------------------
		
	public:
		std::string GetFFName() const;

		/// dipole form factor
		double FF_dipole(double t) const;
	
		/// eletric form factor.
		/// normalized such FF_e(0) = 1, t negative
		double FF_e(double t) const;

		/// magnetic form factor.
		/// normalized such FF_m(0) = 1, t negative
		double FF_m(double t) const;
		
		/// square of the effective form factor.
		/// \param t in GeV^2, negative
		double FF_sq(double t) const;

		/// d/dt of the effective form factor square.
		/// \param t in GeV^2, negative
		double FF_sq_prime(double t) const
		{
			double ep = 1E-5;
			return (FF_sq(t + ep) - FF_sq(t)) / ep;
		}

		//-------------------- interference phases --------------------

		/// Full West-Yennie phase (with alpha factor).
		/// That is the \f$\alpha\Phi\f$ in the decomposition \f$F^{C+H} = F^C e^{i \alpha \Phi} + F^H\f$
		/// Implemented according Eq. (23) in [1]. NB: There is a typo in the formula, there should be
		/// "-" in front of eta in the second term on RHS.
		/// \param t in GeV^2, negative
		TComplex Phi_WY(double t) const;
		
		/// Simplified West-Yennie phase.
		/// Implemented according Eq. (26) in [1].
		TComplex Phi_SWY(double t) const;

		/// Kundrat-Lokajicek phase (with alpha factor).
		/// That is the \f$\alpha\Phi\f$ in the decomposition \f$F^{C+H} = F^C + F^H * e^{i \alpha \Psi}\f$
		/// Implemented according Eq. (26) in [2]. NB: many other publications by Kundrat and Lokajicek have
		/// a wrong sign of the B term!
		/// \param t in GeV^2, negative
		TComplex Psi_KL(double t) const;
		
		/// Kundrat-Lokajicek phase (with alpha factor).
		/// That is the \f$\alpha\Phi\f$ in the decomposition \f$F^{C+H} = F^C + F^H * e^{i \alpha \Psi}\f$
		/// Implemented according Eq. (30) in [3].
		TComplex Psi_Cahn(double t) const;

		/// Interference phase WITH the alpha factor.
		/// returns either \f$-\Phi\f$ or \f$\Psi\f$
		/// \param t in GeV^2, negative
		TComplex Phase(double t) const;

		//-------------------- amplitudes --------------------

		/// pure Coulomb amplitude (PC).
		/// \param t in GeV^2, negative
		TComplex Amp_pure(double t) const;	
		
		TComplex Amp_WY(double t) const;	
		TComplex Amp_SWY(double t) const;	
		TComplex Amp_KL(double t) const;	
		TComplex Amp_Cahn(double t) const;	

		/// total Amplitude according to the choice in `mode'
		TComplex Amp(double t) const;

		//-------------------- standard quantities --------------------

		/// ratio (|KL|^2 - |WY|^2) / |KL|^2.
		/// \param t in GeV^2, negative
		TComplex R(double t) const;
		
		/// for |t| < |cutoff|: (|KL|^2 - |WY|^2) / |KL|^2, otherwise (|KL|^2 - |PH|^2) / |KL|^2.
		/// \param t in GeV^2, negative
		/// \param cutoff in GeV^2, negative
		TComplex R_with_cutoff(double t, double cutoff) const;

		/// ratio (|KL|^2 - |PH|^2 - |PC|^2) / |KL|^2.
		/// \param t in GeV^2, negative
		TComplex Z(double t) const;
		
		/// ratio (|KL|^2 - |PH|^2) / |PH|^2.
		/// \param t in GeV^2, negative
		TComplex C(double t) const;
};

extern CoulombInterference *coulomb;

} // namespace

#endif
