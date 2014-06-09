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

#ifndef _elegent_ppp_model_
#define _elegent_ppp_model_

#include "Model.h"
#include "Math.h"

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
		~PPPModel();
		
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

		bool integ_workspace_initialized;
		unsigned long integ_workspace_size;
		gsl_integration_workspace *integ_workspace;

		static void SetTrajectory(Trajectory &t, double D, double c, double ap, double r2, double s0);

		static TComplex Delta(const Trajectory &, double t);

		/// b in GeV^-1
		virtual TComplex prf0(double b) const;

		static TComplex prf_J0(double b, double *par, const void *obj);
};

} // namespace

#endif
