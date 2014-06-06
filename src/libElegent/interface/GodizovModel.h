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

#ifndef _elegent_godizov_model_
#define _elegent_godizov_model_

#include "Model.h"

#include <gsl/gsl_integration.h>

namespace Elegent
{

/**
 * \brief TODO
 * References:
 *	[1] arXiv:1404.2851v2
 **/
class GodizovModel : public Model
{
	public:
		GodizovModel();
		
		void Configure();

		virtual void Init();

		virtual void Print() const;

		virtual TComplex Amp(double t) const;

		/// b in fm
		virtual TComplex Prf(double b) const;

	protected:
		double De;		// al_P(0) - 1
		double ta_a;
		double Ga_P0;
		double ta_g;
		double s0;

		/// Eq. (2) in [1]
		TComplex delta_t(double t) const;

		static TComplex delta_t_J0(double t, void *vpa);
		static double delta_t_J0_Re(double t, void *vpa);
		static double delta_t_J0_Im(double t, void *vpa);
		
		/// bottom relation from Eq. (1) in [1]
		TComplex delta_b(double b) const;

		/// profile function, b in GeV^-1, see Eq. (1) in [1]
		TComplex prf0(double b) const;

		static TComplex prf_J0(double b, void *vpa);
		static double prf_J0_Re(double b, void *vpa);
		static double prf_J0_Im(double b, void *vpa);

		/// TODO: for GSL
		unsigned long gsl_w_size;
		gsl_integration_workspace *gsl_w;

		/*
		Trajectory pom1, pom2, pom3, oder, regf, rego;
		double s0;
		double precision, upper_bound;


		static void SetTrajectory(Trajectory &t, double D, double c, double ap, double r2, double s0);

		static TComplex Delta(const Trajectory &, double t);

		/// b in GeV^-1
		virtual TComplex prf0(double b) const;
		static TComplex prf_J0(double *b, double *t, const void *obj);
		
		// TODO: for GSL
		static TComplex prf_J0(double b, void *vp);
		static double prf_J0_Re(double b, void *vp);
		static double prf_J0_Im(double b, void *vp);
		*/
};

} // namespace

#endif
