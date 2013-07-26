/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#include <string>
#include <vector>
#include <map>

#include "interface/Constants.h"
#include "interface/AllModels.h"
#include "interface/InterpolationModel.h"

#include "TFile.h"
#include "TGraph.h"

using namespace Elegent;
using namespace std;

//----------------------------------------------------------------------------------------------------

struct AmplitudeGraph
{
	TGraph *re, *im;
};

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: ElegentDistributionSampler [option] [option] ...\n");
	printf("OPTIONS:\n");
	printf("\t-h\t\t\tprint this help and exit\n");
	printf("\t-energy <E>\t\tset incident proton energy (i.e. half of sqrt(s)), in GeV\n");
	printf("\t-pp\t\t\tselect proton-proton interactions\n");
	printf("\t-app\t\t\tselect antiproton-proton interactions\n");
	printf("\t-model-N <number>\tnumber of points for hadronic-model sampling\n");
	printf("\t-model-tmax <value>\tmaximum |t| value (in GeV^2) for hadronic-model sampling\n");
	printf("\t-full-N <number>\tnumber of points in full-range graphs\n");
	printf("\t-full-tmax <value>\tmaximum |t| value (in GeV^2) for full-range graphs\n");
	printf("\t-lowt-N <number>\tnumber of points in low-|t| graphs\n");
	printf("\t-lowt-tmax <value>\tmaximum |t| value (in GeV^2) for low-|t| graphs\n");
	printf("\t-models <string>\tcomma-separated list of model tags\n");
	printf("\t-output <filename>\toutput file name\n");
}

//----------------------------------------------------------------------------------------------------

int InitModels(const string& hadronicModelsString, vector<Model *> &models)
{
	size_t p_curr = 0;
	while (true)
	{
		size_t p_next = hadronicModelsString.find(',', p_curr);
		string tag = (p_next != string::npos) ? hadronicModelsString.substr(p_curr, p_next - p_curr) : hadronicModelsString.substr(p_curr);

		model = NULL;
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
			return 3;
		}

		model->tag = tag;
		models.push_back(model);

		printf("\n>> model with tag `%s':\n", tag.c_str());
		model->Print();

		if (p_next == string::npos)
			break;
		else
			p_curr = p_next + 1;
	}

	return 0;
}

//----------------------------------------------------------------------------------------------------

void SampleModels(const vector<Model *> &models, double t_max, unsigned int N, vector<InterpolationModel *> &models_sampled)
{
	for (unsigned int mi = 0; mi < models.size(); mi++)
	{
		model = models[mi];

		InterpolationModel *ms = new InterpolationModel(N, -t_max, 0.);
		ms->tag = model->tag + ":interpolated";

		for (unsigned int pi = 0; pi < N; pi++)
		{
			double t = ms->GetT(pi);
			ms->SetPoint(pi, model->Amp(t));
		}

		models_sampled.push_back(ms);
	}
}

//----------------------------------------------------------------------------------------------------

void SampleAmplitude(bool logarithmic, unsigned int N, double t_min, double t_max, double t_min_coulomb,
		AmplitudeGraph &a, bool full)
{
	double dt = (t_max - t_min) / (N - 1);
	double xi = pow(t_min/t_max, 1./(N - 1));

	a.re = new TGraph(); a.re->SetName("re");
	a.im = new TGraph(); a.im->SetName("im");
	
	unsigned int idx = 0;
	for (unsigned int pi = 0; pi < N; pi++)
	{
		double t = (logarithmic) ? t_max * pow(xi, pi) : t_max - dt * pi;

		if (!full)
		{
			if (fabs(t) < 1E-6)
				continue;

			if (t < t_min_coulomb)
				continue;
		}

		TComplex v = coulomb->Amp(t);
		a.re->SetPoint(idx, -t, v.Re());
		a.im->SetPoint(idx, -t, v.Im());
		idx++;
	}
}

//----------------------------------------------------------------------------------------------------

void BuildAmplitudes(const vector<InterpolationModel*> &models_sampled,
	bool logarithmic, unsigned int N, double t_min, double t_max,
	const vector<CoulombInterference::ciMode> &amplitudeModes,
	double t_min_coulomb,
	AmplitudeGraph &amplitude_pc, vector< map<CoulombInterference::ciMode, AmplitudeGraph> > &amplitudes)
{
	// PC amplitude
	coulomb->mode = CoulombInterference::mPC;
	SampleAmplitude(logarithmic, N, t_min, t_max, t_min_coulomb, amplitude_pc, false);

	// other amplitudes
	amplitudes.clear();
	amplitudes.resize(models_sampled.size());
	for (unsigned int mi = 0; mi < models_sampled.size(); mi++)
	{
		model = models_sampled[mi];

		map<CoulombInterference::ciMode, AmplitudeGraph> &as = amplitudes[mi];

		for (unsigned int cii = 0; cii < amplitudeModes.size(); cii++)
		{
			coulomb->mode = amplitudeModes[cii];

			AmplitudeGraph &a = as[amplitudeModes[cii]];
			SampleAmplitude(logarithmic, N, t_min, t_max, t_min_coulomb, a, (amplitudeModes[cii] == CoulombInterference::mPH));
		}
	}
}

//----------------------------------------------------------------------------------------------------

void WriteOneAmplitudeGraphs(const AmplitudeGraph &a)
{
	a.re->Write("amplitude_re");
	a.im->Write("amplitude_im");

	// build differential and cumulative cross-sections
	TGraph *phase = new TGraph(); phase->SetName("phase");
	TGraph *rho = new TGraph(); rho->SetName("rho");
	TGraph *B = new TGraph(); B->SetName("B");
	TGraph *dif = new TGraph(); dif->SetName("differential cross-section");
	TGraph *cum = new TGraph(); cum->SetName("cumulative cross-section");

	double *at = a.re->GetX();
	double *ar = a.re->GetY();
	double *ai = a.im->GetY();

	double S=0., prev_t=0., prev_cs=0.;
	for (int i = 0; i < a.re->GetN(); i++)
	{
		TComplex A(ar[i], ai[i]);
		double cs = cnts->sig_fac * A.Rho2();

		dif->SetPoint(i, at[i], cs);
		phase->SetPoint(i, at[i], A.Theta());
		rho->SetPoint(i, at[i], 1./tan(A.Theta()));

		if (i > 0)
			B->SetPoint(i-1, (at[i] + prev_t)/2., (log(prev_cs) - log(cs)) / (at[i] - prev_t));
		
		if (i > 0)
			S += (cs + prev_cs) * (at[i] - prev_t) / 2.;
		
		cum->SetPoint(i, at[i], S);

		prev_cs = cs;
		prev_t = at[i];
	}

	phase->Write();
	rho->Write();
	B->Write();
	dif->Write();
	cum->Write();
}

//----------------------------------------------------------------------------------------------------

void WriteCRZGraphs(const AmplitudeGraph &amp_pc, const AmplitudeGraph &amp_ph,
		const AmplitudeGraph &amp_kl, const AmplitudeGraph &amp_swy)
{
	TGraph *C = new TGraph(); C->SetName("C");
	TGraph *R = new TGraph(); R->SetName("R");
	TGraph *Z = new TGraph(); Z->SetName("Z");

	for (int i = 0; i < amp_kl.re->GetN(); i++)
	{
		double t, pcR, pcI, phR, phI, klR, klI, swyR, swyI;
		amp_kl.re->GetPoint(i, t, klR);
		amp_kl.im->GetPoint(i, t, klI);

		pcR = amp_pc.re->Eval(t);
		pcI = amp_pc.im->Eval(t);
		phR = amp_ph.re->Eval(t);
		phI = amp_ph.im->Eval(t);
		swyR = amp_swy.re->Eval(t);
		swyI = amp_swy.im->Eval(t);

		TComplex ph(phR, phI), pc(pcR, pcI), swy(swyR, swyI), kl(klR, klI);

		double vC = (kl.Rho2() - ph.Rho2()) / ph.Rho2();
		double vZ = (kl.Rho2() - ph.Rho2() - pc.Rho2()) / kl.Rho2();
		double vR = (kl.Rho2() - swy.Rho2()) / kl.Rho2();

		Z->SetPoint(i, t, vZ);
		C->SetPoint(i, t, vC);
		R->SetPoint(i, t, vR);
	}

	C->Write();
	R->Write();
	Z->Write();
}

//----------------------------------------------------------------------------------------------------

void WriteGraphs(const vector<Model *> &models, const AmplitudeGraph &amplitude_pc,
		const vector< map<CoulombInterference::ciMode, AmplitudeGraph> > &amplitudes)
{
	TDirectory *topDir = gDirectory;

	// coulomb graphs
	gDirectory = topDir->mkdir("PC");
	WriteOneAmplitudeGraphs(amplitude_pc);

	for (unsigned int mi = 0; mi < amplitudes.size(); mi++)
	{
		TDirectory *modelDir = topDir->mkdir(models[mi]->tag.c_str());

		for (map<CoulombInterference::ciMode, AmplitudeGraph>::const_iterator mit = amplitudes[mi].begin(); mit != amplitudes[mi].end(); ++mit)
		{
			coulomb->mode = mit->first;
			gDirectory = modelDir->mkdir(coulomb->GetModeString().c_str());

			WriteOneAmplitudeGraphs(mit->second);
		}

		gDirectory = modelDir;
		const map<CoulombInterference::ciMode, AmplitudeGraph> &m = amplitudes[mi];
		map<CoulombInterference::ciMode, AmplitudeGraph>::const_iterator it_ph = m.find(CoulombInterference::mPH);
		map<CoulombInterference::ciMode, AmplitudeGraph>::const_iterator it_kl = m.find(CoulombInterference::mKL);
		map<CoulombInterference::ciMode, AmplitudeGraph>::const_iterator it_swy = m.find(CoulombInterference::mSWY);

		if (it_ph == m.end() || it_kl == m.end() || it_swy == m.end())
		{
			printf("ERROR: some of the PH, KL, SWY amplitudes are missing for model `%s'.\n", models[mi]->tag.c_str());
			continue;
		}

		WriteCRZGraphs(amplitude_pc, it_ph->second, it_kl->second, it_swy->second);
	}
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	// defaults
	double energy = 0.;
	Constants::ParticleMode pMode = Constants::mPP;

	unsigned int model_N = 10000, fullRange_N = 5000, lowt_N = 500;
	double model_t_max = 20., fullRange_t_max = 20., lowt_t_max = 1.;

	string hadronicModelsString = "islam_bfkl,islam_cgc,ppp2,ppp3,bsw,bh,jenkovszky";

	string outputFileName = "";

	// process command line
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-h") == 0)
		{
			PrintUsage();
			return 0;
		}

		if (strcmp(argv[i], "-energy") == 0)
		{
			if (argc-1 > i)
				energy = atof(argv[++i]);
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

		if (strcmp(argv[i], "-model-N") == 0)
		{
			if (argc-1 > i)
				model_N = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-model-tmax") == 0)
		{
			if (argc-1 > i)
				model_t_max = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-full-N") == 0)
		{
			if (argc-1 > i)
				fullRange_N = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-full-tmax") == 0)
		{
			if (argc-1 > i)
				fullRange_t_max = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-lowt-N") == 0)
		{
			if (argc-1 > i)
				lowt_N = atoi(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-lowt-tmax") == 0)
		{
			if (argc-1 > i)
				lowt_t_max = atof(argv[++i]);
			continue;
		}

		if (strcmp(argv[i], "-models") == 0)
		{
			if (argc-1 > i)
				hadronicModelsString = argv[++i];
			continue;
		}


		if (strcmp(argv[i], "-output") == 0)
		{
			if (argc-1 > i)
				outputFileName = argv[++i];
			continue;
		}

		printf("ERROR: unrecognised parameter `%s'.\n", argv[i]);
		PrintUsage();
		return 10;
	}

	// test input
	bool stop = false;
	if (energy == 0.)
	{
		printf("ERROR: energy not specified.\n");
		stop = true;
	}

	if (outputFileName.empty())
	{
		printf("ERROR: output not specified.\n");
		stop = true;
	}

	if (stop)
	{
		PrintUsage();
		return 1;
	}

	// print input
	printf(">> ElegentDistributionSampler > input:\n");
	printf("\tenergy = %E\n", energy);
	printf("\tparticle mode %u\n", pMode);

	printf("\tmodel sampling: samples = %u, t_max = %.2E\n", model_N, model_t_max);
	printf("\tfull-range plots: samples = %u, t_max = %.2E\n", fullRange_N, fullRange_t_max);
	printf("\tlow-|t| plots: samples = %u, t_max = %.2E\n", lowt_N, lowt_t_max);

	printf("\tmodels = %s\n", hadronicModelsString.c_str());
	printf("\toutput = %s\n", outputFileName.c_str());

	// initialise constants etc.
	Constants::Init(2.*energy, pMode);
	cnts->Print();
	coulomb->Print();

	// is model_t_max high enough
	if (model_t_max < coulomb->T)
	{
		printf("ERROR: model_t_max = %.2E is lower than CoulombInterference::T = %.2E\n", model_t_max, coulomb->T);
		return 2;
	}
	double t_min_coulomb = -(model_t_max - coulomb->T);
	
	// inialise models
	vector<Model *> models;
	if (InitModels(hadronicModelsString, models) != 0)
		return 3;

	// sample the models
	vector<InterpolationModel *> models_sampled;
	SampleModels(models, model_t_max, model_N, models_sampled);

	// prepare output
	TFile *outF = new TFile(outputFileName.c_str(), "recreate");

	// a trick to save E, since of->WriteObject(&E, "cmsEnergy") doesn't work
	TGraph *data = new TGraph();
	data->SetName("data");
	data->SetPoint(0, 0., energy);
	data->Write();

	// select amplitude to be generated
	vector<CoulombInterference::ciMode> amplitudeModes;
	amplitudeModes.push_back(CoulombInterference::mPH);
	//amplitudeModes.push_back(CoulombInterference::mWY);
	amplitudeModes.push_back(CoulombInterference::mSWY);
	amplitudeModes.push_back(CoulombInterference::mKL);

	// build full-range amplitudes
	AmplitudeGraph amplitude_pc_full;
	vector< map<CoulombInterference::ciMode, AmplitudeGraph> > amplitudes_full;
	BuildAmplitudes(models_sampled, false, fullRange_N, -fullRange_t_max, 0., amplitudeModes,
		t_min_coulomb, amplitude_pc_full, amplitudes_full);

	// write full-range graphs
	gDirectory = outF->mkdir("full range");
	WriteGraphs(models, amplitude_pc_full, amplitudes_full);

	// build low-|t| amplitudes
	AmplitudeGraph amplitude_pc_lowt;
	vector< map<CoulombInterference::ciMode, AmplitudeGraph> > amplitudes_lowt;
	BuildAmplitudes(models_sampled, true, lowt_N, -lowt_t_max, -1E-5, amplitudeModes,
		t_min_coulomb, amplitude_pc_lowt, amplitudes_lowt);
	
	// build low-|t| graphs
	gDirectory = outF->mkdir("low |t|");
	WriteGraphs(models, amplitude_pc_lowt, amplitudes_lowt);

	delete outF;
	return 0;
}
