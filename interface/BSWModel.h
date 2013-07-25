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
 **/
class BSWModel : public Model 
{
	public:
		/// a Regge trajectory
		struct Trajectory
		{
			double C, b, a, ap, sig;
			void Init(double _C, double _b, double _a, double _ap, double _sig)
				{ C=_C; b=_b; a=_a; ap=_ap; sig = _sig; }
		};

		/// the modes of this model
		enum ModeType { mPomReg, mPom, mReg };

		/// the Pom+Reg mode still broken
		BSWModel() : Model("bsw", "BSW", "Bourrely-Soffer-Wu", mPom) {}
		
		~BSWModel() {}

		/// the Pom+Reg mode still broken
		void Init(ModeType _mode = mPom, bool _presampled = true);

		virtual void Print() const;
		virtual std::string GetModeString() const;

		virtual TComplex Amp(double t) const;	
		virtual TComplex Prf(double b) const;	

	public:
		/// flag whether the presampled mode is on
		bool presampled;

		/// the pomeron exchange parameters
		double c, cp, a, f, m1, m2, asq, m1sq, m2sq;

		/// the 3 Regge trajectories
		Trajectory A2, rho, omega;

		/// temporary constants: Omega0 = S0*F + R0 / s / regge_fac;
		TComplex regge_fac;

		TComplex amp_fac;
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

		static TComplex prf0_J0_presampled(double *b, double *q, const void *obj);

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
