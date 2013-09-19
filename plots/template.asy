import root;
import pad_layout;

//----------------------------------------------------------------------------------------------------

xSizeDef = 8cm;

//----------------------------------------------------------------------------------------------------

string model_tags[];
pen model_pens[];
string model_labels[];

void AddModel(string tag, string label, pen p)
{
	model_tags.push(tag);
	model_labels.push(label);
	model_pens.push(p);
}

AddModel("islam_bfkl", "Islam et al.~(HP)", black);
AddModel("islam_cgc", "Islam et al.~(LxG)", black+dashed);
AddModel("ppp3", "Petrov et al.~(3P)", red);
AddModel("ppp2", "Petrov et al.~(2P)", red+dashed);
AddModel("bsw", "Bourrely et al.", blue);
AddModel("bh", "Block et al.", heavygreen);
AddModel("jenkovszky", "Jenkovszky et al.", magenta);
//AddModel("", " et al.");

//----------------------------------------------------------------------------------------------------

void DrawOptimized(real jump_tol, rObject obj, pen p, string label = "")
{
	int N = obj.iExec("GetN");
	guide g;

	real prev_y = 0;
	for (int i = 0; i < N; ++i)
	{
		real[] xa = {0};
		real[] ya = {0};
		obj.vExec("GetPoint", i, xa, ya);

		real x = xa[0];
		real y = ya[0];

		real de_y = y - prev_y;

		if (i == 0 || abs(de_y) < jump_tol)
		{
			g = g--Scale((x, y));
		} else {
			draw(g, p, label);
			label = "";	// not to multiply labels in legend
			g = Scale((x, y));
		}

		prev_y = y;
	}
	
	draw(g, p, label);
}

//----------------------------------------------------------------------------------------------------

file f_out;

void InitOutput()
{
	f_out = output("include.html");
}

//----------------------------------------------------------------------------------------------------

frame legend_frame;

void PlotAllModels(string input_file, string template, string legend_title,
	bool optimize=false, real opt_jump=0)
{
	for (int mi : model_tags.keys)
	{
		string obj_path = replace(template, "<model>", model_tags[mi]);
		if (optimize)
			DrawOptimized(opt_jump, rGetObj(input_file, obj_path, search=false), model_pens[mi], model_labels[mi]);
		else
			draw(rGetObj(input_file, obj_path, search=false), "l", model_pens[mi], model_labels[mi]);
	}

	legend_frame = BuildLegend(legend_title);
}

//----------------------------------------------------------------------------------------------------

void FinalizeFile(string name, string html_label)
{
	NewPad(false);
	add(legend_frame);

	string sub_dir = "$input_tag";	// $
	
	write(f_out, "        <li><a href=\"" + sub_dir + "/" + name + ".pdf\">" + html_label + "</a></li>", endl);

	GShipout(name);
}

//----------------------------------------------------------------------------------------------------

void GeneratePlots(string mode, string energy, string input_dir, string input_tag)
{
	string mode_html = mode;
	if (mode == "app")
		mode_html = "p&#773;p";

	write(f_out, "<h1>" + mode_html + ", &#8730;s = " + energy + "</h1>", endl);
	
	write(f_out, "    <div style=\"margin-left:3em\">",endl);
	
	write(f_out, "    <h3>t-distributions</h3>",endl);

	// -------------------- t-distributions --------------------

	string input_file = input_dir+"/t-distributions," + input_tag + ".root";
	write(f_out, "    <ul>",endl);

	// ---------- hadronic amplitudes ----------	

	NewPad("$|t|\ung{GeV^2}$", "$\Re F^{\rm H}$");
	scale(Linear, Linear(true));
	// TODO: way to easily add grid?
	PlotAllModels(input_file, "full range/<model>/PH/amplitude_re", "hadronic amplitudes");

	NewPad("$|t|\ung{GeV^2}$", "$\Im F^{\rm H}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "full range/<model>/PH/amplitude_im", "hadronic amplitudes");

	FinalizeFile("t_dist,amp,full_range", "hadronic amplitude, <i>full range</i>");
	
	NewPad("$|t|\ung{GeV^2}$", "$\Re F^{\rm H}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/PH/amplitude_re", "hadronic amplitudes");

	NewPad("$|t|\ung{GeV^2}$", "$\Im F^{\rm H}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/PH/amplitude_im", "hadronic amplitudes");

	FinalizeFile("t_dist,amp,low_t", "hadronic amplitude, <i>low |t|</i>");
	
	// ---------- hadronic differentical cross-section ----------

	NewPad("$|t|\ung{GeV^2}$", "$\d\si / \d t \ung{mb/GeV^2}$");
	scale(Linear, Log(true));
	PlotAllModels(input_file, "full range/<model>/PH/differential cross-section", "hadronic differential cross-section");
	FinalizeFile("t_dist,diff_cs,full_range", "hadronic differential cross-section, <i>full range</i>");
	
	NewPad("$|t|\ung{GeV^2}$", "$\d\si / \d t \ung{mb/GeV^2}$");
	scale(Log, Log(true));
	PlotAllModels(input_file, "low |t|/<model>/PH/differential cross-section", "hadronic differential cross-section");
	FinalizeFile("t_dist,diff_cs,low_t", "hadronic differential cross-section, <i>low |t|</i>");
	
	// ---------- hadronic differential slope ----------
	
	NewPad("$|t|\ung{GeV^2}$", "$B(t) \equiv {\d\over \d t} \log{\d\si\over\d t} \ung{GeV^{-2}}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "full range/<model>/PH/B", "hadronic slope");
	FinalizeFile("t_dist,B,full_range", "hadronic slope, <i>full range</i>");

	NewPad("$|t|\ung{GeV^2}$", "$B(t) \equiv {\d\over \d t} \log{\d\si\over\d t} \ung{GeV^{-2}}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/PH/B", "hadronic slope");
	FinalizeFile("t_dist,B,low_t", "hadronic slope, <i>low |t|</i>");

	// ---------- hadronic phase ----------

	NewPad("$|t|\ung{GeV^2}$", "$\arg F^{\rm H}(t)$");
	scale(Linear, Linear);
	PlotAllModels(input_file, "full range/<model>/PH/phase", "hadronic phase", true, 1);
	ylimits(-3.141593, +3.141593, Crop);
	FinalizeFile("t_dist,phase,full_range", "hadronic phase, <i>full range</i>");

	NewPad("$|t|\ung{GeV^2}$", "$\arg F^{\rm H}(t)$");
	scale(Log, Linear);
	PlotAllModels(input_file, "low |t|/<model>/PH/phase", "hadronic phase", true, 1);
	ylimits(-3.141593, +3.141593, Crop);
	FinalizeFile("t_dist,phase,low_t", "hadronic phase, <i>low |t|</i>");

	// ---------- hadronic rho ----------

	NewPad("$|t|\ung{GeV^2}$", "$\rh(t) \equiv {\Re F^{\rm H}\over \Im F^{\rm H}}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "full range/<model>/PH/rho", "rho parameter", true, 8);
	ylimits(-10, 10, Crop);
	FinalizeFile("t_dist,rho,full_range", "hadronic rho, <i>full range</i>");

	NewPad("$|t|\ung{GeV^2}$", "$\rh(t) \equiv {\Re F^{\rm H}\over \Im F^{\rm H}}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/PH/rho", "rho parameter", true, 8);
	ylimits(-10, 10, Crop);
	FinalizeFile("t_dist,rho,low_t", "hadronic rho, <i>low |t|</i>");

	// TODO: something for the C+H ??
	
	// ---------- interference C ----------

	NewPad("$|t|\ung{GeV^2}$", "$C(t) \equiv {|F^{\rm C+H}|^2 - |F^{\rm H}|^2 \over |F^{\rm H}|^2}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "full range/<model>/C", "influence of Coulomb interaction");
	//ylimits(-10, 10, Crop);
	FinalizeFile("t_dist,C,full_range", "influence of Coulomb interaction, <i>full range</i>");
	
	NewPad("$|t|\ung{GeV^2}$", "$C(t) \equiv {|F^{\rm C+H}|^2 - |F^{\rm H}|^2 \over |F^{\rm H}|^2}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/C", "influence of Coulomb interaction");
	//ylimits(-10, 10, Crop);
	FinalizeFile("t_dist,C,low_t", "influence of Coulomb interaction, <i>low |t|</i>");

	// ---------- interference Z ----------
	
	NewPad("$|t|\ung{GeV^2}$", "$Z(t) \equiv {|F^{\rm C+H}|^2 - |F^{\rm H}|^2 - |F^{\rm C}|^2 \over |F^{\rm C+H}|^2}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "full range/<model>/Z", "importance of the interference term");
	//ylimits(-10, 10, Crop);
	FinalizeFile("t_dist,Z,full_range", "importance of the interference term, <i>full range</i>");

	NewPad("$|t|\ung{GeV^2}$", "$Z(t) \equiv {|F^{\rm C+H}|^2 - |F^{\rm H}|^2 - |F^{\rm C}|^2 \over |F^{\rm C+H}|^2}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/Z", "importance of the interference term");
	//ylimits(-10, 10, Crop);
	FinalizeFile("t_dist,Z,low_t", "importance of the interference term, <i>low |t|</i>");

	// ---------- interference R ----------
	
	NewPad("$|t|\ung{GeV^2}$", "$R(t) \equiv {|F^{\rm KL}|^2 - |F^{\rm WY}|^2 \over |F^{\rm KL}|^2}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "full range/<model>/R", "difference between SWY and KL formulae");
	//ylimits(-10, 10, Crop);
	FinalizeFile("t_dist,R,full_range", "difference between SWY and KL formulae, <i>full range</i>");
	
	NewPad("$|t|\ung{GeV^2}$", "$R(t) \equiv {|F^{\rm KL}|^2 - |F^{\rm WY}|^2 \over |F^{\rm KL}|^2}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/R", "difference between SWY and KL formulae");
	//ylimits(-10, 10, Crop);
	FinalizeFile("t_dist,R,low_t", "difference between SWY and KL formulae, <i>low |t|</i>");
	
	write(f_out, "    </ul>",endl);
	// -------------------- b-distributions --------------------
	
	write(f_out, "    <h3>b-distributions</h3>");
	write(f_out, "    <ul>",endl);
	
	string input_file = input_dir+"/b-distributions," + input_tag + ".root";

	// ---------- hadronic amplitudes ----------	

	NewPad("$b\ung{fm}$", "$\Re T^{\rm H}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "<model>/prf_re", "hadronic amplitudes");

	NewPad("$b\ung{fm}$", "$\Im T^{\rm H}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "<model>/prf_im", "hadronic amplitudes");

	FinalizeFile("b_dist,amp", "hadronic amplitude");
	
	write(f_out, "    </ul>",endl);
	write(f_out, "    </div>",endl);
}

//----------------------------------------------------------------------------------------------------

void FinalizeOutput()
{
}

//----------------------------------------------------------------------------------------------------

InitOutput();
GeneratePlots("$mode", "$energy", "$input_dir", "$input_tag");
FinalizeOutput();
