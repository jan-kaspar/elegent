/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

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
\defgroup Generator Generator
\brief A MC package to generate p-p and p-anti p elastic scattering events.

Generator is ELastic Event GENeraTor. Events are generated according to several phenomenological models, 
see classes IslamModel, PPPModel, BSWModel, BHModel and JenkovszkyModel.

The process is split into two parts
  -# Preparation of cumulative distribution function (CDF). This (lengthy) step is done once and then the results are just reused. 
  The work is done by module CDFBuilder. Some prepared CDFs might be found in data subdirectory.
  -# Generation of events according to a prepared CDF. The work is done by Generator class, which is, for the framework purposes, wrapped
  by the GeneratorSource.
**/


/**
\ingroup Generator
\brief MC generator of proton-proton elastic scattering events.

All internal quantities are in GeV or in mm.
**/
class Generator
{
  public:
    Generator(const std::string &_file, const std::string &_tag, double _t_min, double _t_max, unsigned int _verbosity=1);
    ~Generator() {}

    unsigned int Init();

    static const int PID = 2212;
    static const int ElasticScattering = 91;
    static const int FinalState = 1;
    static const int NullState = 0;

  protected:
    /// cdf source filename
    std::string fileName;

    /// tag of model to be used
    std::string modelTag;

    /// |t| values in GeV^2, bounds for cdf
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
