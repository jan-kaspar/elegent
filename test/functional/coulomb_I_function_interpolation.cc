#include "interface/Constants.h"
#include "interface/CoulombInterference.h"

#include "TFile.h"
#include "TProfile2D.h"

#include <string>
#include <chrono>
#include <iostream>

using namespace std;
using namespace Elegent;

//----------------------------------------------------------------------------------------------------

void TimingTest()
{
	for (double mt = 0; mt <= 20.; mt += 0.02)
	{
		for (double mtp = 0; mtp <= 20.; mtp += 0.02)
		{
			const double t_min_diff = 0.1;	// GeV^2
			if (fabs(mt - mtp) < t_min_diff)
				continue;

			coulomb->I_function(-mt, -mtp);
		}
	}
}

//----------------------------------------------------------------------------------------------------

void RunAccuracyTestSquare(unsigned int n_t_bins, double t_max, const string &name, const string &message)
{
	// timing variables
	chrono::time_point<chrono::system_clock> t_start, t_end;
	chrono::duration<double> t_duration;

	// run loop
	TProfile2D *p_rel_eff_tp_vs_t = new TProfile2D(name.c_str(), ";t;t'", n_t_bins, 0., t_max, n_t_bins, 0., t_max);

	cout << "* starting accuracy test: " << message << endl;
	t_start = chrono::system_clock::now();
	for (unsigned int i = 0; i < n_t_bins; i++)
	{
		for (unsigned int j = 0; j < n_t_bins; j++)
		{
			const double w = t_max / n_t_bins;
			const double mt = w * (0.5 + i);
			const double mtp = w * (0.5 + j);

			const double t_min_diff = 1E-5;	// GeV^2
			if (fabs(mt - mtp) < t_min_diff)
				continue;

			const double I0 = -2.*M_PI / fabs(mt - mtp);
	
			coulomb->useIFunctionInterpolator = false;
			const double I_exa = coulomb->I_function(-mt, -mtp);

			coulomb->useIFunctionInterpolator = true;
			const double I_int = coulomb->I_function(-mt, -mtp);

			const double rI_exa = I_exa / I0;
			const double rI_int = I_int / I0;

			const double rel_diff = (rI_int - rI_exa) / rI_exa;

			//printf("%.2E\n", rel_diff);
			p_rel_eff_tp_vs_t->Fill(mt, mtp, fabs(rel_diff));
		}
	}
	t_end = chrono::system_clock::now();

	t_duration = t_end - t_start;
	cout << "    done in " << t_duration.count() << " s" << endl;

	p_rel_eff_tp_vs_t->Write();
}

//----------------------------------------------------------------------------------------------------

void RunAccuracyTestDiagonal(unsigned int n_m_bins, double m_max,
		unsigned int n_d_bins, double d_min, double d_max,
		const string &name, const string &message)
{
	// timing variables
	chrono::time_point<chrono::system_clock> t_start, t_end;
	chrono::duration<double> t_duration;

	// run loop
	TProfile2D *p_rel_eff_d_vs_m = new TProfile2D(name.c_str(), ";m = (|t'| + |t|)/2;d = |t'| - |t|", n_m_bins, 0., m_max, n_d_bins, d_min, d_max);

	cout << "* starting accuracy test: " << message << endl;
	t_start = chrono::system_clock::now();
	for (unsigned int i = 0; i < n_m_bins; i++)
	{
		for (unsigned int j = 0; j < n_d_bins; j++)
		{
			const double m_w = m_max / n_m_bins;
			const double m = m_w * (0.5 + i);

			const double d_w = (d_max - d_min) / n_d_bins;
			const double d = d_min + d_w * (0.5 + j);

			if (fabs(d/2) > m)
				continue;

			const double mt = m - d/2;
			const double mtp = m + d/2;

			const double mt_max = 20.;
			if (mt > mt_max || mtp > mt_max)
				continue;

			const double t_min_diff = 1E-5;	// GeV^2
			if (fabs(mt - mtp) < t_min_diff)
				continue;

			const double I0 = -2.*M_PI / fabs(mt - mtp);

			coulomb->useIFunctionInterpolator = false;
			const double I_exa = coulomb->I_function(-mt, -mtp);

			coulomb->useIFunctionInterpolator = true;
			const double I_int = coulomb->I_function(-mt, -mtp);

			const double rI_exa = I_exa / I0;
			const double rI_int = I_int / I0;

			const double rel_diff = (rI_int - rI_exa) / rI_exa;

			//printf("%.2E\n", rel_diff);
			p_rel_eff_d_vs_m->Fill(m, d, fabs(rel_diff));
		}
	}
	t_end = chrono::system_clock::now();

	t_duration = t_end - t_start;
	cout << "    done in " << t_duration.count() << " s" << endl;

	p_rel_eff_d_vs_m->Write();
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// init Elegent
	Constants::Init(2*6500, cnts->mPP);
    cnts->Print();

	coulomb->ffType = coulomb->ffPuckett;
	coulomb->precision = 1E-5;
	coulomb->Print();

	// prepare output
	TFile *f_out = TFile::Open("coulomb_I_function_interpolation.root", "recreate");

	// timing variables
	chrono::time_point<chrono::system_clock> t_start, t_end;
	chrono::duration<double> t_duration;

	// prepare interpolator
	cout << "* initiating interpolator" << endl;
	t_start = chrono::system_clock::now();
	coulomb->InitIFunctionInterpolator(20., 300);
	t_end = chrono::system_clock::now();

	t_duration = t_end - t_start;
	cout << "    done in " << t_duration.count() << " s" << endl;

	// test exact evaluation
	cout << "* starting timing loop with exact evaluation" << endl;
	t_start = chrono::system_clock::now();
	coulomb->useIFunctionInterpolator = false;
	TimingTest();
	t_end = chrono::system_clock::now();

	t_duration = t_end - t_start;
	cout << "    done in " << t_duration.count() << " s" << endl;

	// test interpolator
	cout << "* starting timing loop with interpolator" << endl;
	t_start = chrono::system_clock::now();
	coulomb->useIFunctionInterpolator = true;
	TimingTest();
	t_end = chrono::system_clock::now();

	t_duration = t_end - t_start;
	cout << "    done in " << t_duration.count() << " s" << endl;

	// accuracy tests
	RunAccuracyTestSquare(100, 20., "p_rel_eff_tp_vs_t.full", "square, full");
	RunAccuracyTestSquare(100, 0.3, "p_rel_eff_tp_vs_t.low_t", "square, low t");

	RunAccuracyTestDiagonal(100, 20., 100, -0.1, +0.1, "p_rel_eff_d_vs_m.full", "diagonal, full");
	RunAccuracyTestDiagonal(100, 0.3, 100, -0.1, +0.1, "p_rel_eff_d_vs_m.low_t", "diagonal, low_t");

	RunAccuracyTestDiagonal(100, 20., 100, 0.5, 0.6, "p_rel_eff_d_vs_m.offdiagonal", "off-diagonal");

	// clean up
	delete f_out;

	return 0;
}
