// HEADER

#include "interface/IslamModel.h"
#include "interface/Math.h"
#include "interface/Constants.h"

using namespace std;
using namespace Elegent;

//#define DEBUG 1

//----------------------------------------------------------------------------------------------------

IslamModel::IslamModel() : Model("islam", "islam(uninitialized)", "Islam et al."),
  precision(1E-11), precision_t(1E-3), upper_bound(50.), upper_bound_t(-15.)
{
}

//-------------------------------------- DIFFRACTION AMPLITUDE ------------------------------------

TComplex IslamModel::GammaD(double b) const
{
  /// profile function Gamma_D^+
  /// b... impact parameter in  fm
  return 1. / (1. + TComplex::Exp((b - R) / a)) + 1. / (1. + TComplex::Exp((-b - R) / a)) - 1.;
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::GammaD_J0(double *b, double *t, const void *obj)
{
  /// b[0] ... impact parameter in fm
  /// t[0] ... t in GeV^2
  return ((IslamModel *)obj)->GammaD(b[0]) * b[0] * TMath::BesselJ0(b[0]*sqrt(-t[0]));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_diff(double t) const
{
#ifdef DEBUG
  printf(">> IslamModel::T_diff\n");
#endif

  /// t ... t
  return Diff_fac * CmplxInt(this, GammaD_J0, 0, upper_bound, &t, precision);
}

//----------------------------------------- CORE AMPLITUDE ----------------------------------------

double IslamModel::F_sq(double t)  const
{
  /// formfactor, t ... t
  return beta * sqrt(M_sq - t) * TMath::BesselK1(beta * sqrt(M_sq - t));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_core(double t) const
{
#ifdef DEBUG
  printf(">> IslamModel::T_core\n");
#endif

  /// t ... t
  return Core_fac * Hard_fac * F_sq(t) / (M_sq - t);
}


//---------------------------------------- QUARK AMPLITUDE ----------------------------------------

double IslamModel::I_integral(double qt, double al) const
{
  double ap = qt/2./al;
  double a = sqrt(ap * ap + 1.);
  double ival;
  
  // for small qt values just substitute limit value 16/3
  if (qt > 1E-10)
    ival = (  2./a/a + 1./ap/ap - 3.*ap*ap/a/a/a/a  ) / a/ap * log(a + ap)  -  1./a/a/ap/ap + 3./a/a/a/a;
  else ival = 16./3.;

  return 1./8. /al/al/al/al * ival;
}

//----------------------------------------------------------------------------------------------------

double IslamModel::F_cal_integ(double *x, double *par, const void *obj)
{
  double &qt = par[0];
  double &n = par[1];
  double &omega = par[2];
  double &m0sq = par[3];
  double al_sq = m0sq/4. + cnts->M_sq * x[0] * x[0];
  double al = sqrt(al_sq);
  return exp((1. + n * omega) * log(x[0])) / al_sq * ((IslamModel *)obj)->I_integral(qt, al);
}

//----------------------------------------------------------------------------------------------------

double IslamModel::F_cal(int n, double qt, double omega, double m0sq) const
{
  double par[4];
  par[0] = qt;
  par[1] = n;
  par[2] = omega;
  par[3] = m0sq;
  return cnts->M * exp(2.5 * log(m0sq)) / 8. / cnts->pi * DoubleInt(this, F_cal_integ, 0., 1., par, 1E-9);
}

//----------------------------------------------------------------------------------------------------

double IslamModel::T_qq_integ(double *bArr, double *par, const void *obj)
{
  double &b = bArr[0];
  double &q = par[0]; // par[0] ... q
  double &n = par[1]; // par[1] ... n
  return b * TMath::BesselJ0(b * q) * pow( TMath::BesselK0(b / ((IslamModel *)obj)->r0) , n);
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_qq(int n, double t) const
{
  double q = sqrt(fabs(t));

  if (n == 1)
    return Quark_fac / (-t + 1./r0/r0);
  
  if (n == 2) {
    if (q < 1E-2)
      return i * Quark_fac*Quark_fac * r0*r0*r0/4.;  // limit
    else
      return i * Quark_fac*Quark_fac * 2. * asinh(q * r0 / 2.) / q / sqrt( fabs(t) + 4./r0/r0 );
  }
  
  // general formula
  double par[2];
  par[0] = q;
  par[1] = n;
  // correct -i
  return -i  / 2. / TMath::Factorial(n) * TComplex::Power(Quark_const, n) * DoubleInt(this, T_qq_integ, 0., 30., par, 1E-9);
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_quark(double t) const
{
#ifdef DEBUG
  printf(">> IslamModel::T_quark\n");
#endif

  /// t ... t
  double qt = sqrt(-t * (-t / cnts->t_min + 1.));

  TComplex sum = 0.;
  for (int j = 1; j <= qqMaxOrder; j++) {
    double F = F_cal(j, qt, omega, m0sq);
    sum += F*F * T_qq(j, t);
  }

  return Hard_fac * sum;
}

//------------------------------------ CGC AMPLITUDE -------------------------------------------------

double IslamModel::T_cgc_integ(double *bArr, double *par, const void *obj)
{
  double &b = bArr[0];
  double &q = par[0]; // par[0] ... q
  double &n = par[1]; // par[1] ... n
  double &m_c = ((IslamModel *)obj)->m_c;
  return b * TMath::BesselJ0(b * q) * pow(exp(-b * m_c) * m_c*m_c * (1. + b*m_c) / 3., n);
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_cgc_n(int n, double t) const
{
  double q = sqrt(fabs(t));

  if (n == 1)
    return i * cgc_fac / pow(1. - t/m_c/m_c, 2.5);
  
  // general formula
  double par[2];
  par[0] = q;
  par[1] = n;
  return i  * pow(-2., n - 1) / TMath::Factorial(n) * TComplex::Power(cgc_fac, n) * DoubleInt(this, T_cgc_integ, 0., 30., par, 1E-9);
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_cgc(double t) const
{
#ifdef DEBUG
  printf(">> IslamModel::T_cgc\n");
#endif

  /// t ... t
  double qt = sqrt(-t * (-t / cnts->t_min + 1.));

  TComplex sum = 0.;
  for (int j = 1; j <= cgcMaxOrder; j++) {
    double F = F_cal(j, qt, lambda, m0sq);
    sum += F*F * T_cgc_n(j, t);
  }

  return Hard_fac * sum;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------- FULL AMPLITUDE ----------------------------------------

TComplex IslamModel::Amp(double t) const
{  
#ifdef DEBUG
  printf(">> IslamModel::amp, mode = %i\n", mode);
#endif

  /// t ... t

  switch (mode) {
    case mFullQuark:  return T_diff(t) + T_core(t) + T_quark(t);
    case mFullCGC:    return T_diff(t) + T_core(t) + T_cgc(t);
    case mDiff:       return T_diff(t);
    case mCore:       return T_core(t);
    case mQuark:      return T_quark(t);
    case mCGC:        return T_cgc(t);
    case mDiffCore:   return T_diff(t) + T_core(t);
    default:
      printf("!!! mode = -1 !!!\n\n"); return 0.;
  }
}


//----------------------------------------------------------------------------------------------------
//---------------------------------------- PROFILE FUNCTIONS --------------------------------------

TComplex IslamModel::T_core_J0(double *t, double *b, const void *obj)
{
  /// b[0] ... impact parameter in fm
  /// t[0] ... t in GeV^2
  return ((IslamModel *)obj)->T_core(t[0]) * TMath::BesselJ0(b[0]*sqrt(-t[0]));
}

//----------------------------------------------------------------------------------------------------
TComplex IslamModel::T_diff_J0(double *t, double *b, const void *obj)
{
  /// b[0] ... impact parameter in fm
  /// t[0] ... t in GeV^2
  return ((IslamModel *)obj)->T_diff(t[0]) * TMath::BesselJ0(b[0]*sqrt(-t[0]));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::T_quark_J0(double *t, double *b, const void *obj)
{
  /// b[0] ... impact parameter in fm
  /// t[0] ... t in GeV^2
  return ((IslamModel *)obj)->T_quark(t[0]) * TMath::BesselJ0(b[0]*sqrt(-t[0]));
}

//----------------------------------------------------------------------------------------------------

TComplex IslamModel::Prf(double b) const
{
  // conversion from fm to GeV^-1
  /*
  b /= cnts->hbarc; 
  if (mode == 0) return Diff_fac_profile * GammaD(b) + (  CmplxInt(this, T_core_J0, upper_bound_t, 0, &b, precision_t) + CmplxInt(this, T_quark_J0, upper_bound_t, 0, &b, 1E-9)  ) / 4. / cnts->p_cms / cnts->sqrt_s ;
  if (mode == 1) return Diff_fac_profile * GammaD(b);
  if (mode == 2) return CmplxInt(this, T_core_J0, upper_bound_t, 0, &b, precision_t) / 4. / cnts->p_cms / cnts->sqrt_s;
  if (mode == 3) return CmplxInt(this, T_quark_J0, upper_bound_t, 0, &b, precision_t) / 4. / cnts->p_cms / cnts->sqrt_s;
  if (mode == 4) return Diff_fac_profile * GammaD(b) + CmplxInt(this, T_core_J0, upper_bound_t, 0, &b, precision_t) / 4. / cnts->p_cms / cnts->sqrt_s;
  */
  return 0;
}


//----------------------------------------------------------------------------------------------------
//------------------------------------------- INITIALISATION --------------------------------------

TComplex CEF(double a, double b, double c)
{
// returns a + b / (s exp(-i pi/2))^c = a + b / (-i s)^c = a + b (-i s)^(-c)
  return a + b * TComplex::Power(-i * cnts->s, -c);
}

//----------------------------------------------------------------------------------------------------

void IslamModel::InitBase(double R0, double R1, double a0, double a1, double be, double _m)
{
  R = TComplex(R0 + R1*cnts->ln_s, -R1*cnts->pi/2);
  a = TComplex(a0 + a1*cnts->ln_s, -a1*cnts->pi/2);

  beta = be; M_sq = _m*_m;
}

//----------------------------------------------------------------------------------------------------

void IslamModel::InitStage1(double pi_g_mod, double g_arg, double pi_hga, double hth, double GaMD_Re, double GaMD_Im)
{
  // enum g(s) first
  Diff_fac = pi_g_mod / cnts->pi * TComplex::Exp(i * g_arg);
  
  if (cnts->pMode == cnts->mPP) {
    Hard_fac = cnts->s / cnts->pi * pi_hga * TComplex::Exp(i * hth) * (1 + GaMD_Re +i * GaMD_Im - (1. - TComplex::Exp(-R/a)) / (1. + TComplex::Exp(-R/a)) * Diff_fac); 
    Core_fac = -1;
  }
  
  if (cnts->pMode == cnts->mAPP) {
    Hard_fac = + cnts->s / cnts->pi * pi_hga * TComplex::Exp(i * hth) * (1 - GaMD_Re -i * GaMD_Im - (1. - TComplex::Exp(-R/a)) / (1. + TComplex::Exp(-R/a)) * Diff_fac); 
    Core_fac = +1;
  }
  
  // Difffraction coeficient f_D, in Diff_fac is stored g(s)
  Diff_fac_profile = i / 2. * Diff_fac;
  Diff_fac = i * cnts->p_cms * cnts->sqrt_s * Diff_fac;
  
  name = "Islam(ST1)";
}

//----------------------------------------------------------------------------------------------------

void IslamModel::InitStage2(double et0, double c0, double si, double la0, double d0, double al, double hga0, double hga1, double hsi)
{
  Diff_fac_profile = (1. - CEF(et0, c0, si)) * (1. + TComplex::Exp(-R/a)) / (1. - TComplex::Exp(-R/a));
  Diff_fac = i * cnts->sqrt_s * cnts->p_cms * Diff_fac_profile;
  Diff_fac_profile = Diff_fac_profile * i/2.;

  if (cnts->pMode == cnts->mPP) {
    Hard_fac = cnts->s * CEF(hga0, hga1, hsi) * ( CEF(et0, c0, si) + i*CEF(la0, -d0, al) );
    Core_fac = -1;
  }

  if (cnts->pMode == cnts->mAPP) {
    Hard_fac = cnts->s * CEF(hga0, hga1, hsi) * ( CEF(et0, c0, si) - i*CEF(la0, -d0, al) );
    Core_fac = +1;
  }

  name = "Islam(ST2)";
}

//----------------------------------------------------------------------------------------------------

void IslamModel::InitQQ(double _tgaqq, double _omega, double _r0, double _m0sq)
{
  m0sq = _m0sq; r0 = _r0; omega = _omega;
  Quark_fac = i * _tgaqq * TComplex::Power(-i * cnts->s, _omega);
  Quark_const = -2. * _tgaqq * TComplex::Power(-i * cnts->s, _omega);
  qqMaxOrder = 1;
  name = "Islam(ST3)";
}

//----------------------------------------------------------------------------------------------------

void IslamModel::InitCGC(double _tgagg, double _lambda, double _m_c, double _m0sq)
{
  m0sq = _m0sq; lambda = _lambda; m_c = _m_c;
  cgc_fac = _tgagg * TComplex::Power(-i * cnts->s, lambda);
  cgcMaxOrder = 1;
  name = "Islam(ST4)";
}

//----------------------------------------------------------------------------------------------------

void IslamModel::DirectInit(double Rr, double Ri, double ar, double ai, double Dr, double Di, double Hr, double Hi)
{
  // STAGE 2
  R = TComplex(Rr, Ri);
  a = TComplex(ar, ai);
  Diff_fac = TComplex(Dr, Di);
  Hard_fac = TComplex(Hr, Hi);
  name = "Islam(ST2direct)";
}

//----------------------------------------------------------------------------------------------------
//------------------------------------------- PRINTING --------------------------------------------

void IslamModel::Print() const
{
  printf(">> IslamModel::Print\n");
  printf("old\n");
  double v1 = R.Im(); v1 = -v1 * 2 / cnts->pi;
  double v0 = R.Re(); v0 -= v1 * cnts->ln_s;
  printf("\tR0=%f\n\tR1=%f", v0, v1);
  v1 = a.Im(); v1 = -v1 * 2 / cnts->pi;
  v0 = a.Re(); v0 -= v1 * cnts->ln_s;
  printf("\n\ta0=%f\n\ta1=%f", v0, v1);
  printf("\n\tRe f_D=%f\n\tIm f_D=%f\n\tRe f_H=%f\n\tIm f_H=%f\n", Diff_fac.Re(), Diff_fac.Im(), Hard_fac.Re(), Hard_fac.Im());
  
  printf("diffraction variables\n");
  printf("\tR: Re=%E, Im=%E\n", R.Re(), R.Im());
  printf("\ta: Re=%E, Im=%E\n", a.Re(), a.Im());
  printf("\tDiff_fac_profile: Re=%E, Im=%E\n", Diff_fac_profile.Re(), Diff_fac_profile.Im());
  printf("\tDiff_fac: Re=%E, Im=%E\n", Diff_fac.Re(), Diff_fac.Im());
  printf("\tHard_fac: Re=%E, Im=%E\n", Hard_fac.Re(), Hard_fac.Im());

  printf("core scattering variables\n");
  printf("\tbeta = %E\n", beta);
  printf("\tM_sq = %E\n", M_sq);
  printf("\tCore_fac = %E\n", Core_fac);

  printf("quark-quard scattering variables\n");
  printf("\tm0sq = %E\n", m0sq);
  printf("\tr0 = %E\n", r0);
  printf("\tomega = %E\n", omega);
  printf("\tQuark_fac: Re=%E, Im=%E\n", Quark_fac.Re(), Quark_fac.Im());
  printf("\tQuark_const: Re=%E, Im=%E\n", Quark_const.Re(), Quark_const.Im());
  printf("\tqqMaxOrder = %i\n", qqMaxOrder);
  printf("\tlambda = %E\n", lambda);
  printf("\tm_c = %E\n", m_c);
  printf("\tcgc_fac: Re=%E, Im=%E\n", cgc_fac.Re(), cgc_fac.Im());
  printf("\tcgcMaxOrder = %i\n", cgcMaxOrder);

  printf("integration variables\n");
  printf("\tprecision = %E\n", precision);
  printf("\tprecision_t = %E\n", precision_t);
  printf("\tupper_bound = %E\n", upper_bound);
  printf("\tupper_bound_t = %E\n", upper_bound_t);
}

//----------------------------------------------------------------------------------------------------

string IslamModel::GetModeString() const
{
  switch (mode) {
    case mDiff:       return "diffraction";
    case mCore:       return "core";
    case mQuark:      return "quark";
    case mCGC:        return "CGC";
    case mDiffCore:   return "core+diffraction";
    case mFullQuark:  return "core+diffraction+quark";
    case mFullCGC:    return "core+diffraction+CGC";
    default:          return "unknown";
  }
}

// vim: ft=cmscpp
