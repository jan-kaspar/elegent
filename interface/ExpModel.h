/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_exp_model_
#define _elegent_exp_model_

#include "Model.h"

#include <cstring>

namespace Elegent
{


/**
 * \ingroup Elegent
 * \brief Exponential model of p-p and p-anti p elastic scattering.
 * For reference pourposes only.
 **/
class ExpModel : public Model
{
  public:
    ExpModel();
    
    virtual void Print() const;

    virtual std::string GetModeString() const
      { return "basic"; }
  
    virtual TComplex Amp(double t) const;
    virtual TComplex Prf(double b) const;

  public:
    double a, b1, b2, p0, p1;
};

} // namespace

#endif
