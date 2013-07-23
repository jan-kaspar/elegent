// HEADER

#include "interface/InterpolationModel.h"

#include "TGraph.h"

namespace Elegent
{

//#define DEBUG 1

//----------------------------------------------------------------------------------------------------

InterpolationModel::InterpolationModel() :
  Model("interpolation", "", ""),
  t_min(0.), t_max(0.), re(new TGraph()), im(new TGraph())
{
}

//----------------------------------------------------------------------------------------------------

InterpolationModel::~InterpolationModel()
{
  delete re;
  delete im;
}

//----------------------------------------------------------------------------------------------------

void InterpolationModel::Print() const
{
  printf(">> InterpolationModel::Print\n");
  printf("\tpoints stored %u\n", re->GetN());

  for (int idx = 0; idx < re->GetN(); idx++) {
    double t, r, i;
    re->GetPoint(idx, t, r);
    im->GetPoint(idx, t, i);

    printf("\t\t%i\t%E\t%E\t%E\n", idx, t, r, i);
  }
}

//----------------------------------------------------------------------------------------------------

TComplex InterpolationModel::Prf(double) const
{
  printf(">> InterpolationModel::prf > Not implemented.\n");
  return 0;
}

//----------------------------------------------------------------------------------------------------

TComplex InterpolationModel::Amp(double t) const
{
  /// convert |t| to t
  double mt = -t;

  // to disable extrapolation
  if (mt < t_min || mt > t_max)
    return TComplex(0, 0);

  return TComplex(re->Eval(mt, NULL, ""), im->Eval(mt, NULL, ""));
}

//----------------------------------------------------------------------------------------------------

void InterpolationModel::AddPoint(double t, double r, double i)
{
  /// convert t to |t|
  double mt = -t;

  if (re->GetN() == 0)
  {
    t_min = mt;
    t_max = mt;
  } else {
    if (mt < t_min)
      t_min = mt;
    if (mt > t_max)
      t_max = mt;
  }

  re->SetPoint(re->GetN(), mt, r);
  im->SetPoint(im->GetN(), mt, i);
}    

} // namespace
