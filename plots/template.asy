import root;
import pad_layout;

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

file f_out;

void InitOutput()
{
	f_out = output("include.html");
}

//----------------------------------------------------------------------------------------------------

frame legend_frame;

void PlotAllModels(string input_file, string template, string legend_title)
{
	for (int mi : model_tags.keys)
	{
		string obj_path = replace(template, "<model>", model_tags[mi]);
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
	write(f_out, "<h2>&#8730;s = " + energy + "</h2>", endl);
	
	write(f_out, "    <h3>t-distributions</h3>",endl);

	// -------------------- t-distributions --------------------

	string input_file = input_dir+"/t-distributions," + input_tag + ".root";

	// ---------- hadronic amplitudes, full range ----------	

	NewPad("$|t|\ung{GeV^2}$", "$\Re F^{\rm H}$");
	scale(Linear, Linear(true));
	// TODO: way to easily add grid?
	PlotAllModels(input_file, "full range/<model>/PH/amplitude_re", "hadronic amplitudes");

	NewPad("$|t|\ung{GeV^2}$", "$\Im F^{\rm H}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "full range/<model>/PH/amplitude_im", "hadronic amplitudes");

	FinalizeFile("t_dist,amp,full_range", "hadronic amplitude, full range");
	
	// ---------- hadronic amplitudes, low |t| ----------	

	NewPad("$|t|\ung{GeV^2}$", "$\Re F^{\rm H}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/PH/amplitude_re", "hadronic amplitudes");

	NewPad("$|t|\ung{GeV^2}$", "$\Im F^{\rm H}$");
	scale(Log, Linear(true));
	PlotAllModels(input_file, "low |t|/<model>/PH/amplitude_im", "hadronic amplitudes");

	FinalizeFile("t_dist,amp,low_t", "hadronic amplitude, low |t|");
	
	// ---------- differentical cross-section, hadronic, full range ----------

	NewPad("$|t|\ung{GeV^2}$", "$\d\si / \d t \ung{mb/GeV}$");
	scale(Linear, Log(true));
	PlotAllModels(input_file, "full range/<model>/PH/differential cross-section", "hadronic differential cross-section");

	FinalizeFile("t_dist,diff_cs,full_range", "hadronic differential cross-section, full range");

	// ---------- phase, hadronic, full range ----------
	// ---------- rho, hadronic, full range ----------
	// ---------- differential slope, hadronic, full range ----------

	// TODO: something for the C+H ??
	
	// ---------- interference C, full range ----------
	// ---------- interference R, full range ----------
	// ---------- interference Z, full range ----------
	
	// -------------------- b-distributions --------------------
	
	write(f_out, "    <h3>b-distributions</h3>");
	
	string input_file = input_dir+"/b-distributions," + input_tag + ".root";

	// ---------- hadronic amplitudes ----------	

	// TODO: check the units of b !! Also add to ROOT file graph!
	NewPad("$b\ung{fm}$", "$\Re T^{\rm H}$");
	scale(Linear, Linear(true));
	PlotAllModels(input_file, "<model>/prf_re", "hadronic amplitudes");

	FinalizeFile("b_dist,amp", "hadronic amplitude");
}

//----------------------------------------------------------------------------------------------------

void FinalizeOutput()
{
}

//----------------------------------------------------------------------------------------------------

InitOutput();
GeneratePlots("$mode", "$energy", "$input_dir", "$input_tag");
FinalizeOutput();
