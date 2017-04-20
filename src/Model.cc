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

#include "interface/Model.h"

using namespace std;

namespace Elegent
{

string Model::CompileShortLabel() const
{
	string out = shortLabel.name;
 
  	if (!shortLabel.mode.empty() && shortLabel.mode.compare("full"))
		out += ":" + shortLabel.mode;

  	if (!shortLabel.variant.empty())
		out += " (" + shortLabel.variant + ")";

	out += " [" + shortLabel.version + "]";

	return out;
}

//----------------------------------------------------------------------------------------------------

string Model::CompileFullLabel() const
{
	string out = fullLabel.name;
 
  	if (!fullLabel.mode.empty() && fullLabel.mode.compare("full"))
		out += ":" + fullLabel.mode;

  	if (!fullLabel.variant.empty())
		out += " (" + fullLabel.variant + ")";

	out += " [" + fullLabel.version + "]";

	return out;
}

//----------------------------------------------------------------------------------------------------

Model *model = NULL;

}
