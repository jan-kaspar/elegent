/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_model_
#define _elegent_model_

#include "TComplex.h"

#include <string>

namespace Elegent
{

/**
 * \brief The base class for hadronic models of (anti)proton-proton elastic scattering.
 **/
class Model
{
	public:
		std::string tag;	///< short model name
		std::string name;	///< model name, version, etc.
		std::string label;	///< long model name

		signed char mode;	///< model mode

		Model(const std::string& _tag="", const std::string &_name="", const std::string &_label="", signed char _mode=-1) :
			tag(_tag), name(_name), label(_label), mode(_mode) {}

		virtual ~Model() {}

		/// prints model info
		virtual void Print() const	=0;

		/// returns model mode as a string
		virtual std::string GetModeString() const =0;

		///\brief amplitude, t in GeV^-2, t < 0
		/// Normalisation is such that
		///   dsigma/dt = (\hbar c)^2 * \pi / (s p^2) * |Amp(t)|^2
		/// Differential cross-section can be obtained as
		///   Constants::sig_fac * |Amp(t)|^2
		virtual TComplex Amp(double t) const =0;

		///\brief profile function, b in fm
		/// Normalisation is such that
		///   Amp(t) = 2 p \sqrt{s} \int db b Prf() J_0(b \sqrt{-t})
		virtual TComplex Prf(double b) const =0;
};

/// actual model
extern Model *model;

} // namespace

#endif
