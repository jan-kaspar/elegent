/********************************************************************************

    Copyright 2013 Jan Kašpar

    This file is part of Elegent (http://elegent.hepforge.org/).

    Elegent is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Elegent is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Elegent.  If not, see <http://www.gnu.org/licenses/>.
 
********************************************************************************/

#include "interface/ModelFactory.h"

using namespace std;

namespace Elegent
{

ModelFactory::ModelFactory()
{
	BHModel *bh = new BHModel();
	bh->Configure();
	model_map[bh->CompileShortLabel()] = bh;

	BSWModel *bsw = new BSWModel();
	bsw->Configure(BSWModel::mPomReg);
	model_map[bsw->CompileShortLabel()] = bsw;

	DLModel *dl = new DLModel();
	dl->Configure();
	model_map[dl->CompileShortLabel()] = dl;

	FerreiraModel *ferreira = new FerreiraModel();
	ferreira->Configure();
	model_map[ferreira->CompileShortLabel()] = ferreira;
	
	GodizovModel *godizov = new GodizovModel();
	godizov->Configure();
	model_map[godizov->CompileShortLabel()] = godizov;
	
	IslamModel *islam_lxg = new IslamModel();
	islam_lxg->Configure(IslamModel::vLxG, IslamModel::mFull);
	model_map[islam_lxg->CompileShortLabel()] = islam_lxg;
	
	IslamModel *islam_hp = new IslamModel();
	islam_hp->Configure(IslamModel::vHP, IslamModel::mFull);
	model_map[islam_hp->CompileShortLabel()] = islam_hp;

	IslamModel2015 *islam15 = new IslamModel2015();
	islam15->Configure(IslamModel2015::mFull);
	model_map[islam15->CompileShortLabel()] = islam15;
	
	JenkovszkyModel *jenkovszky = new JenkovszkyModel();
	jenkovszky->Configure();
	model_map[jenkovszky->CompileShortLabel()] = jenkovszky;

	PPPModel *ppp2 = new PPPModel();
	ppp2->Configure(PPPModel::v2P);
	model_map[ppp2->CompileShortLabel()] = ppp2;
	
	PPPModel *ppp3 = new PPPModel();
	ppp3->Configure(PPPModel::v3P);
	model_map[ppp3->CompileShortLabel()] = ppp3;
}

//----------------------------------------------------------------------------------------------------

void ModelFactory::PrintList() const
{
	printf(">> ModelFactory::PrintList > available models:\n");

	for (map<std::string, Model*>::const_iterator it = model_map.begin(); it != model_map.end(); ++it)
	{
		printf("\t%s : %s\n", it->first.c_str(), it->second->CompileFullLabel().c_str());
	}
}

//----------------------------------------------------------------------------------------------------

Model* ModelFactory::MakeInstance(const std::string &tag, bool callInit) const
{
	// look for tag in model map
	map<string, Model*>::const_iterator it = model_map.find(tag);
	Model* model = (it == model_map.end()) ? NULL : it->second;
	
	// if not found print all possibilities
	if (model == NULL)
	{
		printf("ERROR in ModelFactory::MakeInstance: model tag `%s' not available\n", tag.c_str());
		PrintList();
		return NULL;
	}

	// initialise model
	if (callInit)
		model->Init();

	return model;
}

} // namespace
