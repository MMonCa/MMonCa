/*
 * Author: ignacio.martin@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "AnnealCmd.h"
#include "domains/Global.h"
#include "io/ParameterManager.h"
#include "kernel/Domain.h"
#include "lkmc/LatticeParam.h"
#include "lkmc/EpiGasParam.h"

#include <iomanip>

using std::string;
using std::map;

namespace IO {

AnnealCmd::AnnealCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv)
{
}
	
int AnnealCmd::operator()()
{
	float annealTime = getFloat("time");
	float temp       = getFloat("temp");
	bool bDepth      = specified("depth");
	float events     = (specified("events")? getFloat("events") : 0);
	float depth      = bDepth? getFloat("depth") : 0;
	std::map<std::string, float> epitaxy = (specified("epitaxy") ? getFloatMap("epitaxy") :
			std::map<std::string, float>());

	//setting epitaxial gases in the proper places
	for(unsigned dom = 0; dom < Domains::global()->getDomains(); ++dom)
		for(Kernel::M_TYPE mt = 0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
			if(Domains::global()->PM()->getMaterial(mt)._bEpitaxy)
				Domains::global()->getDomain(dom)->_pEGPar[mt]->processLine(_pTcl, Domains::global()->getDomain(dom), mt, temp + 273.15, epitaxy);

	Domains::global()->setTempK(temp + 273.15); //recomputes all the rates

	LOWMSG("Annealing the sample for " << annealTime << " seconds at " << temp+273.15 <<
			"K (" << temp << "ºC)");

	Domains::global()->anneal(annealTime, bDepth, depth, events);
	Domains::global()->meshReport();
	Domains::global()->defectReport();
	Domains::global()->eventReport();
	Domains::global()->reactionReport();

	return TCL_OK;	
}

}

