/**************************************************
 * This file is a part of the Elegent package:
 * 	http://elegent.hepforge.org/
 *************************************************/

#ifndef _elegent_modelfactory_
#define _elegent_modelfactory_

#include "IslamModel.h"
#include "PPPModel.h"
#include "BSWModel.h"
#include "BHModel.h"
#include "JenkovszkyModel.h"
#include "ExpModel.h"

#include <string>

namespace Elegent
{

/**
 * \brief TODO
 **/
class ModelFactory
{
	public:
		static Model* MakeStandardInstance(const std::string &tag, bool prf_presampled = true);
};


} // namespace

#endif
