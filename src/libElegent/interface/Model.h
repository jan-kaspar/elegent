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
		/// collection of strings that describe a model instance
		struct Label
		{
			std::string name, variant, version, mode;
		};

		/// full label (e.g. for figure legend)
		Label fullLabel;	

		/// short label (e.g. for object names in ROOT files)
		Label shortLabel;

		/// compiles a human readable string from fullLabel
		std::string CompileFullLabel() const;

		/// compiles a human readable string from shortLabel
		std::string CompileShortLabel() const;

		Model() {}

		virtual ~Model() {}

		/// sets up model parameters and data members
		virtual void Init() =0;

		/// prints model info
		virtual void Print() const	=0;

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

		///\brief sets the presampling option, if available
		/// This option determines whether calling Init would presample relevant distributions
		/// (typically b distributions) for faster Amp calls. This behaviour is, in contrary,
		/// undesirable for evaluating s distributions.
		virtual void ForcePresampling(bool /*value*/) {}
};

/// current model
extern Model *model;

} // namespace

#endif
