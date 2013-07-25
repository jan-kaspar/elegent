/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_jenkovszky_model_
#define _elegent_jenkovszky_model_

#include "Model.h"

namespace Elegent
{

/**
 * \brief The model of Jenkovszky et al.
 * According to L. Jenkovszky, A. Lengyel, D. Lontkovskyi:
 * The Pomeron and Odderon in elastic, inelastic and total cross sections at the LHC.
 * arXiv:1105.1202 (2011)
 **/
class JenkovszkyModel : public Model
{
	public:
		JenkovszkyModel() : Model("jenkovszky", "jenkovszky (default)", "Jenkovszky et al.")
			{ Init(); }

		void Init();

		virtual void Print() const;
		virtual std::string GetModeString() const
			{ return "default"; }

		virtual TComplex Amp(double t) const;
		virtual TComplex Prf(double b) const;

	protected:
		/// pomeron parameters
		double a_P, b_P, de_P, al1_P, ep_P, s_P;

		/// odderon parameters
		double a_O, b_O, de_O, al1_O, s_O;

		/// omega parameters
		double a_om, b_om, s_om, al0_om, al1_om;

		/// f parameters
		double a_f, b_f, s_f, al0_f, al1_f;
};

} // namespace

#endif
