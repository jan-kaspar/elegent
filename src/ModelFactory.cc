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
		IslamBFKL->Init(IslamModel::mFullQuark);
		IslamBFKL->name = "Islam (HP)";
		model = IslamBFKL;
	}

	if (tag.compare("islam_cgc") == 0)
	{
		IslamModel *IslamCGC = new IslamModel();
		IslamCGC->Init(IslamModel::mFullCGC);
		IslamCGC->name = "Islam (LxG)";
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
