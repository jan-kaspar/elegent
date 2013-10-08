/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_bsw_model_
#define _elegent_bsw_model_

#include "Model.h"

#include <vector>

//#define DEBUG

namespace Elegent
{

/**
 * \brief Bourelly, Soffer and Wu model of p-p and p-anti p elastic scattering.
 * References:
 *	[1] BOURRELY C., SOFFER, J. and WU, T. T., Phys. Rev. D19 (1979) 3249
 *	[2] BOURRELY C., SOFFER, J. and WU, T. T., Nucl. Phys. B247 (1984) 15
 *	[3] BOURRELY C., SOFFER, J. and WU, T. T., Eur. Phys. J. C28 (2003) 97-105
 *	[4] BOURRELY C., SOFFER, J. and WU, T. T., Eur. Phys. J. C71 (2011) 1601
 **/
class BSWModel : public Model 
{
	public:
		/// a Regge trajectory
		struct Trajectory
		{
			/// the parameters from Eq. (7) in [3]
			double C, b, a, ap, signature;

			/// Non-documented overall sign of the amplitude. These signs make
			/// the difference between pp and app reactions (see the Init method).
			double sign;

			void Init(double _C, double _b, double _a, double _ap, double _signature)
			{
				C=_C; b=_b; a=_a; ap=_ap; signature = _signature;
				sign = 0.;
			}
		};

		/// available modes
		enum ModeType
		{
			mPomReg,	///< both Pomeron and Reggeon contributions
			mPom,		///< only Pomeron contribution
			mReg		///< only Reggeon contribution
		} mode;

		BSWModel();
		
		~BSWModel() {}

		void Configure(ModeType _mode = mPomReg, bool _presampled = true);
		
		virtual void Init();

		virtual void Print() const;

		virtual TComplex Amp(double t) const;	
		virtual TComplex Prf(double b) const;	

	public:
		/// flag whether the presampled mode is on
		bool presampled;

		/// accuracy level (use true for differential cross-section, false is sufficient for total cross-section)
		bool highAccuracy;

		/// the pomeron exchange parameters
		double c, cp, a, f, m1, m2, asq, m1sq, m2sq;

		/// the 3 Regge trajectories
		Trajectory A2, rho, omega;

		/// constants to resolve ambiguities in the source papers
		TComplex regge_fac;		///< Omega0 = S0*F + R0 / s / regge_fac
		signed int k_u;			///< u = -|u| exp(i * (2 k pi - pi))
		signed int k_lnu;		///< ln u = |ln u| exp(i * (al + 2 k_lnu pi)), al = atan2(Im ln u, Re ln u) in (-pi, +pi)

		double upper_bound_t, precision_t;
		double upper_bound_b, precision_b;

		/// \tilde F(t), Eq. (4)
		double Ft(double t) const;

		/// generic Regge trajectory amplitude, Eq. (7)
		TComplex Rt(Trajectory tr, double t) const;

		/// the sum of allowed Regge trajctories (A2, rho, omega)
		TComplex R0t(double t) const;

		/// S_0(s), Eq(3)
		TComplex S0(double t) const;
		
		/// S_0(0)
		TComplex S00;

		/// the Bessel transform of \Omega_0(s, b) in Eq. (2)
		TComplex Omega0t(double t) const;

		static TComplex Omega0t_J0(double *t, double *b, const void *obj);

		/// \Omega_0(s, b)
		TComplex Omega0b(double b) const;

		/// the profile function with b in GeV^-1
		TComplex prf0(double b) const;

		static TComplex prf0_J0(double *b, double *q, const void *obj);

		/// the sampling-step size
		double data_db;

		/// the number of sampled points
		unsigned int data_N;

		/// the sampled real and imaginary values of prf0(b) 
		std::vector<double> data_re, data_im;

#ifdef DEBUG
		std::vector<double> data_b;
#endif

		/// samples the prf0 function
		void BuildSample(unsigned int samples);

		/// interpolates (linearly) the sample at point b
		TComplex SampleEval(double b);
};

} // namespace

#endif
