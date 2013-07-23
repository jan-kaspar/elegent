// HEADER

#ifndef _elegent_interpolation_model_
#define _elegent_interpolation_model_

#include "Model.h"

class TGraph;

namespace Elegent
{

/**
 * \ingroup Elegent
 * \brief Model that interpolates stored amplitude points.
 **/
class InterpolationModel : public Model
{
  protected:
    /// the lower boundary of the supported interval
    double t_min;

    /// the upper boundary of the supported interval
    double t_max;

    /// graph |t| versus the real part of the amplitude
    TGraph *re; //->

    /// graph |t| versus the imaginary part of the amplitude
    TGraph *im; //->

  public:
    InterpolationModel();
    virtual ~InterpolationModel();

    virtual std::string GetModeString() const
      { return "basic"; }

    virtual void Print() const;             ///< prints model info

    virtual TComplex Amp(double t) const;   ///< amplitude, t in GeV^-2, t < 0
    virtual TComplex Prf(double b) const;   ///< profile function, b in fm

    /// to store a new point
    /// IMPORTANT: points MUST be inserted in ascending |t| order
    void AddPoint(double t, double r, double i);    
};

} // namespace

#endif
