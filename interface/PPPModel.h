/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_ppp_model_
#define _elegent_ppp_model_

#include "Model.h"

namespace Elegent
{

/**
 * \brief Predazzi, Petrov and Prokudin model of p-p and p-anti p elastic scattering.
 * References:
 *	[1] PETROV, V. A. and PROKUDIN, A. V., Eur. Phys. J. C23 (2002) 135–143
 **/
class PPPModel : public Model
{
	public:
		struct Trajectory
		{
			double D, c, ap, r2, rho2;
			TComplex gamma;
		};

		/// available variants
		enum VariantType
		{
			v2P,	///< with 2 Pomerons
			v3P		///< with 3 Pomerons
		} variant;
		
		PPPModel();
		
		void Configure(VariantType v);

		virtual void Init();

		virtual void Print() const;

		virtual TComplex Amp(double t) const;

		/// b in fm
		virtual TComplex Prf(double b) const;

	protected:
		Trajectory pom1, pom2, pom3, oder, regf, rego;
		double s0;
		double precision, upper_bound;

		static void SetTrajectory(Trajectory &t, double D, double c, double ap, double r2, double s0);

		static TComplex Delta(const Trajectory &, double t);

		/// b in GeV^-1
		virtual TComplex prf0(double b) const;
		static	TComplex prf_J0(double *b, double *t, const void *obj);
};

} // namespace

#endif
