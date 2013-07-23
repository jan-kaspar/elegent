// HEADER

//#include <cstdio>
//#include <cstdlib>
#include <string>

#include "interface/Constants.h"
#include "interface/AllModels.h"

using namespace Elegent;
using namespace std;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	// TODO
	printf("USAGE: ElegentDistributionSampler\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	// defaults
	double energy = 0.;
	unsigned int samples = 4000;
	double t_max = 20.;
	string outputFileName = "";
	string hadronicModelsString = "islam_bfkl,islam_cgc,ppp2,ppp3,bsw,bh,jenkovszky";
	Constants::ParticleMode pMode = Constants::mPP;

	// process command line
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-energy") == 0)
		{
			if (argc-1 > i)
				energy = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-samples") == 0)
		{
			if (argc-1 > i)
				samples = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-tmax") == 0)
		{
			if (argc-1 > i)
				t_max = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-output") == 0)
		{
			if (argc-1 > i)
				outputFileName = argv[++i];
			continue;
		}

		if (strcmp(argv[i], "-pp") == 0)
		{
			pMode = Constants::mPP;
			continue;
		}

		if (strcmp(argv[i], "-app") == 0)
		{
			pMode = Constants::mAPP;
			continue;
		}
	}

	// test input
	bool stop = false;
	if (energy == 0.)
	{
		printf("ERROR: energy not specified.\n");
		PrintUsage();
		stop = true;
	}

	if (outputFileName.empty())
	{
		printf("ERROR: output not specified.\n");
		PrintUsage();
		stop = true;
	}

	if (stop)
		return 1;

	// print input
	printf("particle mode %u\n", pMode);
	printf("energy = %E\n", energy);
	printf("t_max = %E\n", t_max);
	printf("samples = %u\n", samples);
	printf("output = %s\n", outputFileName.c_str());
	printf("models = %s\n", hadronicModelsString.c_str());

	// initialise constants etc.
	Constants::Init(2.*energy, pMode);
	cnts->Print();
	coulomb->Print();

	// initilise the selected models
	size_t p_curr = 0;
	while (true)
	{
		size_t p_next = hadronicModelsString.find(',', p_curr);
		string tag = (p_next != string::npos) ? hadronicModelsString.substr(p_curr, p_next - p_curr) : hadronicModelsString.substr(p_curr);

		model = NULL;
		if (tag.compare("islam_bfkl") == 0)
		{
			IslamModel *IslamBFKL = new IslamModel();
			IslamBFKL->InitBase(  2.77,   0.0491, 0.245,   0.126,   3.075,   0.801);
			IslamBFKL->InitStage2(  0.0844,  0.0,  2.7,  0.727,  13.,  0.246,  1.53,  0.,    1.46);
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
			IslamCGC->InitBase(  2.77,   0.0491, 0.245,   0.126,   3.075,   0.801);
			IslamCGC->InitStage2(  0.0844,  0.0,  2.7,  0.727,  13.,  0.246,  1.53,  0.,    1.46);
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
			BSW->Init();
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
		{
			printf("ERROR: model tag `%s' not recognised\n", tag.c_str());
			return 2;
		}

		model->tag = tag;

		printf("\n>> model with tag `%s':\n", tag.c_str());
		model->Print();

		if (p_next == string::npos)
			break;
		else
			p_curr = p_next + 1;
	}

	return 0;
}
