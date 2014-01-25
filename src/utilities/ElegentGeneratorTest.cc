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

#include <cstdio>
#include <cstdlib>

#include "interface/Generator.h"

#include "HepMC/GenEvent.h"

using namespace Elegent;
using namespace std;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: ElegentGeneratorTest <file name> <model path> <t_min> <t_max> <number of events>\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
	if (argc != 6)
	{
		PrintUsage();
		return 1;
	}

	Generator generator(argv[1], argv[2], atof(argv[3]), atof(argv[4]));	
	if (generator.Init() != 0)
		return 2;

	unsigned int N = atoi(argv[5]);
	for (unsigned int i = 0; i < N; i++)
	{
		GenEvent* gEv = new GenEvent();
		gEv->set_event_number(i);
		generator.Generate(gEv);
		gEv->print();
		delete gEv;
	}

	return 0;
}
