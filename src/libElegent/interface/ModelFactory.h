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

#ifndef _elegent_modelfactory_
#define _elegent_modelfactory_

#include "BHModel.h"
#include "BSWModel.h"
#include "DLModel.h"
#include "ExpModel.h"
#include "FerreiraModel.h"
#include "GodizovModel.h"
#include "IslamModel.h"
#include "PPPModel.h"
#include "JenkovszkyModel.h"

#include <string>
#include <map>

namespace Elegent
{

/**
 * \brief A class to give list of available models and to create an instance of a model specified by tag.
 **/
class ModelFactory
{
	protected:
		/// map: tag --> model instance
		std::map<std::string, Model*> model_map;

	public:
		ModelFactory();

		void PrintList() const;

		Model* MakeInstance(const std::string &tag, bool callInit = true) const;
};


} // namespace

#endif
