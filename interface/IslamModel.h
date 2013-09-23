/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_islam_model_
#define _elegent_islam_model_

#include "Model.h"

namespace Elegent
{


/**
 * \brief Islam model of p-p and p-anti p elastic scattering.
 * References:
 *	[1] ISLAM, M. M., FEARNLEY, T. and GUILLAUD, J. P., Nuovo Cim. A81 (1984) 737
 *	[2] ISLAM, M. M., INNOCENTE V., FEARNLEY T. and SANGUIETTI, G., Europhys. Lett. 4 (1987) 189-196
 *	[3] ISLAM, M. M., LUDDY, R. J. and PROKUDIN, A. V., Phys. Lett. B605 (2005) 115-122
 *	[4] ISLAM, M. M., LUDDY, R. J. and PROKUDIN, A. V., Int. J. Mod. Phys. A21 (2006) 1-42
 *	[5] ISLAM, M. M., KASPAR, J. and LUDDY, R. J., Mod. Phys. Lett. A24 (2009) 485-496
 **/
class IslamModel : public Model
{
	public:
		/// modes of Islam model
		enum {mDiff, mCore, mQuark, mCGC, mDiffCore, mFullQuark, mFullCGC};

		IslamModel();

		/// initialization methods
		void Init(unsigned int _mode);
		void InitBase(double R0, double R1, double a0, double a1, double be, double _m);
		void InitStage1(double pi_g_mod, double g_arg, double pi_hga, double hth, double GaMD_Re, double GaMD_Im);
		void InitStage2(double et0, double c0, double si, double la0, double d0, double al, double hga0, double hga1, double hsi);
		void InitQQ(double _tgaqq, double _omega, double _r0, double _m0sq);
		void InitCGC(double _tgagg, double _lambda, double _m_c, double _m0sq);

		virtual void Print() const;
		virtual std::string GetModeString() const;

		virtual TComplex Amp(double t) const;
		virtual TComplex Prf(double b) const;

		void SetUnitarizationOrders(int qq, int cgc)
			{ qqMaxOrder = qq; cgcMaxOrder = cgc; }
		
	protected:
		/// diffraction variables
		TComplex R, a, Diff_fac_profile, Diff_fac;
		
		/// hard scattering variables
		TComplex Hard_fac;
	
		/// core scattering variables
		double beta, M_sq, Core_fac;
	
		/// quark-quard scattering variables
		double m0sq, r0, omega;
		TComplex Quark_fac;
		TComplex Quark_const;
		int qqMaxOrder;
			
		/// CGC scattering
		double lambda, m_c;
		TComplex cgc_fac;
		int cgcMaxOrder;
	
		/// integration variables
		double precision, precision_t, upper_bound, upper_bound_t;

		/// amplitude methods
		TComplex GammaD(double b) const;
		static TComplex GammaD_J0(double *b, double *t, const void *obj);
		TComplex T_diff(double t) const;
		
		double F_sq(double t) const;
		TComplex T_core(double t) const;

		double I_integral(double qt, double al) const;
		static double F_cal_integ(double *x, double *qt, const void *obj);
		double F_cal(int n, double qt, double om, double m0sq) const;

		static double T_qq_integ(double *, double *, const void *);
		TComplex T_qq(int n, double t) const;
		TComplex T_quark(double t) const;

		static double T_cgc_integ(double *, double *, const void *);
		TComplex T_cgc_n(int n, double t) const;
		TComplex T_cgc(double t) const;

		/// profile funcion methods
		static TComplex Amp_J0(double *t, double *b, const void *obj);
};

} // namespace

#endif
