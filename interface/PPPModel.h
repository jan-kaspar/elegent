// HEADER

#ifndef _elegent_ppp_model_
#define _elegent_ppp_model_

#include "Model.h"

namespace Elegent
{

/**
 * \ingroup Elegent
 * \brief Predazzi, Petrov and Prokudin model of p-p and p-anti p elastic scattering.
 **/
class PPPModel : public Model
{
	public:
		struct Trajectory {
			double D, c, ap, r2, rho2;
			TComplex gamma;
		};

		/// modes of PPP model
		enum { m2P, m3P };
		
		PPPModel() : Model("ppp", "PPP (uninitialized)", "", -1) {}
		
		void Init(unsigned char m);

		virtual void Print() const;
		virtual std::string GetModeString() const;

		virtual TComplex Amp(double t) const;
		virtual TComplex Prf(double b) const;		 ///< b in fm

	protected:
		Trajectory pom1, pom2, pom3, oder, regf, rego;
		double s0;
		double precision, upper_bound;

		static void SetTrajectory(Trajectory &t, double D, double c, double ap, double r2, double s0);

		TComplex tr_eik(Trajectory, double s, double t) const;


		virtual TComplex prf0(double b) const;		///< b in GeV^-1
		static	TComplex prf_J0(double *b, double *t, const void *obj);
};

} // namespace

#endif
