/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#include "interface/ModelFactory.h"

namespace Elegent
{

Model* ModelFactory::MakeStandardInstance(const std::string &tag, bool prf_presampled)
{
	Model* model = NULL;
	
	if (tag.compare("islam_bfkl") == 0)
	{
		IslamModel *IslamBFKL = new IslamModel();
		IslamBFKL->InitBase(	2.77,	 0.0491, 0.245,	 0.126,	 3.075,	 0.801);
		IslamBFKL->InitStage2(	0.0844,	0.0,	2.7,	0.727,	13.,	0.246,	1.53,	0.,		1.46);
		IslamBFKL->InitQQ(0.03, 0.15, 2., 12.);
		IslamBFKL->InitCGC(0.0056, 0.29, 1.67, 12.);
		IslamBFKL->SetUnitarizationOrders(1, 1);
		IslamBFKL->mode = IslamModel::mFullQuark;
		IslamBFKL->name = "Islam (BFKL)";
		model = IslamBFKL;
	}

	if (tag.compare("islam_cgc") == 0)
	{
		IslamModel *IslamCGC = new IslamModel();
		IslamCGC->InitBase(	2.77,	 0.0491, 0.245,	 0.126,	 3.075,	 0.801);
		IslamCGC->InitStage2(	0.0844,	0.0,	2.7,	0.727,	13.,	0.246,	1.53,	0.,		1.46);
		IslamCGC->InitQQ(0.03, 0.15, 2., 12.);
		IslamCGC->InitCGC(0.0056, 0.29, 1.67, 12.);
		IslamCGC->SetUnitarizationOrders(1, 1);
		IslamCGC->mode = IslamModel::mFullCGC;
		IslamCGC->name = "Islam (CGC)";
		model = IslamCGC;
	}

	if (tag.compare("ppp2") == 0)
	{
		PPPModel *PPP2 = new PPPModel();
		PPP2->Init(PPPModel::m2P);
		model = PPP2;
	}

	if (tag.compare("ppp3") == 0)
	{
		PPPModel *PPP3 = new PPPModel();
		PPP3->Init(PPPModel::m3P);
		model = PPP3;
	}
	
	if (tag.compare("bsw") == 0)
	{
		BSWModel *BSW = new BSWModel();
		BSW->Init(BSWModel::mPom, prf_presampled);
		model = BSW;
	}

	if (tag.compare("bh") == 0)
	{
		BHModel *BH = new BHModel();
		BH->Init();
		model = BH;
	}

	if (tag.compare("jenkovszky") == 0)
	{	
		JenkovszkyModel *Jenkovszky = new JenkovszkyModel();
		Jenkovszky->Init();
		model = Jenkovszky;
	}

	if (model == NULL)
		printf("ERROR in ModelFactory::MakeStandardInstance: model tag `%s' not recognised\n", tag.c_str());

	return model;
}

} // namespace
