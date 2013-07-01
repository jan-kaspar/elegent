// HEADER

#ifndef _elegent_constants_
#define _elegent_constants_

#include "Rtypes.h" // TODO: needed?

#include "TComplex.h"

namespace Elegent
{


/**
 *\ingroup Elegent
 *\brief Set of constants used in Elegent calculations.
 **/
struct Constants
{
  // physics constants
  static double alpha;					///< fine structure constant
  static double proton_mass;			///< GeV
  static double neutron_mass;			///< GeV
  static double hbarc;					///< GeV * fm
  static double sq_hbarc;				///< GeV^2 * mbarn
  										///< sigma/mbarn = sigma/GeV^-2 * sq_hbarc
  static double M;                      ///< abbreviation for proton mass in GeV
  static double M_sq;                   ///< proton mass squared, GeV^2
  static double kappa;                  ///> the anomalous magnetic moment of proton
  
  // mathematics constants
  static double pi;						///< pi
  static double gamma;					///< Euler's constant
  
  // physics data and preenumerations
  double sqrt_s;
  double s;
  double ln_s;
  double p_cms;
  double sig_fac;
  double t_min;
  
  signed char pMode;					    ///< particle mode
  enum {mPP, mAPP};
 
  /// \param W = sqrt(s)
  Constants(double W = 0., char mode = mPP);

  /// \param W = sqrt(s)
  static void Init(double W, char mode);

  /// print the actual values
  void Print();
};


extern TComplex i;
extern Constants *cnts;

} // namespace

#endif

// vim: ft=cmscpp
