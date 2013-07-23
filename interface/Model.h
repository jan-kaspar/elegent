// HEADER

#ifndef _elegent_model_
#define _elegent_model_

#include "TComplex.h"

#include <string>

namespace Elegent
{

/**
 * \ingroup Elegent
 * \brief The base class for models of p-p and p-anti p elastic scattering.
 **/
class Model
{
  public:
    std::string tag;                              ///< short model name
    std::string name;                             ///< model name, version, etc.
    std::string label;                            ///< long model name

    signed char mode;                             ///< model mode

    Model() : tag(""), name(""), label(""), mode(-1) {}

    Model(const std::string& _tag="", const std::string &_name="", const std::string &_label="", signed char _mode=-1) :
      tag(_tag), name(_name), label(_label), mode(_mode) {}

    virtual ~Model() {}

    virtual void Print() const  =0;               ///< prints model info
    virtual std::string GetModeString() const =0; ///< returns model mode as a string

    virtual TComplex Amp(double t) const =0;      ///< amplitude, t in GeV^-2, t < 0
    virtual TComplex Prf(double b) const =0;      ///< profile function, b in fm
};

/// actual model
extern Model *model;

} // namespace

#endif
