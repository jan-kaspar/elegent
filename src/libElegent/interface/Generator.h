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

#ifndef _elegent_elegent_
#define _elegent_elegent_

#include <string>

class TGraph;

namespace HepMC {
	class GenEvent;
}

namespace Elegent
{

/**
 * MC generator of proton-proton elastic scattering events.
 *
 * All internal quantities are in GeV or in mm.
 **/
class Generator
{
	public:
		Generator(const std::string &_file, const std::string &_path, double _t_min, double _t_max, unsigned int _verbosity=1);
		~Generator() {}

		unsigned int Init();

		static const int PID = 2212;
		static const int ElasticScattering = 91;
		static const int FinalState = 1;
		static const int NullState = 0;

	protected:
		/// name of file containing the cumulative distribution function (CDF)
		std::string fileName;

		/// path of the (CDF) within the file
		std::string modelPath;

		/// |t| values in GeV^2, bounds for CDF
		double t_min, t_max;
		
		/// verbosity level (0 = no, 1 = normal, 2 = debug)
		unsigned int verbosity;

		/// [GeV] cms (one) proton energy
		double E_cms;

		/// [GeV] cms proton momentum
		double p_cms;

		/// graph with inverse c.d.f.
		TGraph *icdf;
	
	public:
		/// generates one event provided two random numbers with uniform distribution on (0, 1)
		void GenerateBase(double rn1, double rn2, HepMC::GenEvent* gE);

		/// generates one event, using ROOT random number generator TRandom2
		void Generate(HepMC::GenEvent* gE);
};

} // namespace

#endif
