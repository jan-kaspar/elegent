// HEADER

#ifndef _elegent_bh_model_
#define _elegent_bh_model_

#include "Model.h"

namespace Elegent
{

/**
 * \ingroup Elegent
 * \brief Block-Halzen model of p-p and p-anti p elastic scattering.
 * REFERENCE PAPER: Block et al., Phys. Rev. D 60, 054024 (1999)
 **/
class BHModel : public Model {
  public:
    BHModel() : Model("bh", "Blk-Hlz", "Block-Halzen", -1) {}

    void Init();

    virtual std::string GetModeString() const
      { return "basic"; }

    virtual void Print() const;

    virtual TComplex Amp(double t) const;
    virtual TComplex Prf(double b) const;

  protected:
    /// common parameters
    double m0, s0, al_s, Sigma_gg; 

    /// parameters for sigma_gg
    double Cp_gg, epsilon, Ng; 
    double a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5;
    
    /// parameters for sigma_gg
    double C_qg_log; 
    
    /// parameters for sigma_qg
    double C, C_even_regge; 
    
    /// parameters for sigma_qq
    double C_odd; 

    /// effective areas
    double mu_gg, mu_qg, mu_qq, mu_odd;

    /// the integral cross-sections
    TComplex sigma_gg, sigma_qq, sigma_qg, sigma_odd;
    
    /// integration parameters
    double precision, upper_bound;

    double W(double b, double mi) const;
    double sumR1(double s) const;
    double sumR2(double s) const;
    double sumR(double s) const;
    double sumI1(double s) const;
    double sumI2(double s) const;
    double sumI(double s) const;

    /// the full eikonal: Eq. (1) without the leading factor i and Eq. (12)
    TComplex chi_without_i(double b) const;

    TComplex prf0(double b) const;
    static TComplex prf0_J0(double *y, double *t, const void *obj);
};

} // namespace

#endif

// vim: ft=cmscpp
