int check_graph_points(const char *file, const char *obj, unsigned N_exp)
{
	TFile *f = TFile::Open(file);
	printf("file pointer: %p\n", f);
	if (!f)
		return 1;

	TGraph *g = (TGraph *) f->Get(obj);
	printf("graph pointer: %p\n", g);
	if (!g)
		return 2;

	int N = g->GetN();
	printf("graph points: %i\n", N);

	if (N != N_exp)
		return 3;

	return 0;
}
