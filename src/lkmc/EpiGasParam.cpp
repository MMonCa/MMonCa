/*
 * EpiGasParam.cpp
 *
 *  Created on: Oct 22, 2014
 *  Author: ignacio.martin@imdea.org
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

#include <sstream>
#include "EpiGasParam.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "io/FileParameters.h"
#include "lkmc/LatticeParam.h"

using std::map;
using std::string;
using std::stringstream;

namespace LKMC {

EpiGasParam::EpiGasParam()
{
	initialize();
}

void EpiGasParam::initialize()
{
	for(Kernel::P_TYPE pt = 0; pt < Kernel::MAX_IMPURITIES; ++pt)
		_partialPress[pt] = 0;
	for(unsigned i=0; i<MAX_EPI_GASES; ++i)
		_toPType[i] = Kernel::UNDEFINED_TYPE;
}

void EpiGasParam::processLine(Tcl_Interp *pTcl, const Kernel::Domain *pDom, Kernel::M_TYPE mt, double T, const map<string, float> &line)
{
	if(line.size() > MAX_EPI_GASES)
		ERRORMSG("Epitaxy gas: Maximum number of allowed gases is " << MAX_EPI_GASES <<
				" but " << line.size() << " were specified.");
	initialize();
	unsigned index = 0;
	stringstream args;
	for(map<string, float>::const_iterator it=line.begin(); it!=line.end(); ++it)
	{
		Kernel::P_TYPE pt = Domains::global()->PM()->getElementNumber(it->first);
		args << it->first << " " << it->second << " ";
		if(pt == Kernel::UNDEFINED_TYPE)
		{
			WARNINGMSG("Particle " << it->first << " not recognized as an epitaxy gas");
			continue;
		}
		pt = Domains::global()->PM()->getFamily(pt);
		_partialPress[pt] = it->second;
		_toPType[index++] = pt;
	}
	//add needed arguments, but they are null ones!
	for(unsigned i=line.size(); i < MAX_EPI_GASES; ++i)
		args << "none 0.0 ";
	args << T;
	if(pDom->_pLaPar[mt]->_type == LatticeParam::DIAMOND || pDom->_pLaPar[mt]->_type == LatticeParam::DIAMOND2)
		for(unsigned i=0; i < index; ++i)
		{
			std::stringstream ss;
			Kernel::P_TYPE myPt = _toPType[i];
			ss << args.str() << " " << Domains::global()->PM()->getParticleName(mt, myPt);
			_prefPrecursor[myPt] = Domains::global()->getFileParameters()->getFloatProcArgs(
					pTcl,
					Domains::global()->PM()->getMaterialName(mt) + "/Epitaxy/prefactor.precursor",
					ss.str());
		}
}

} /* namespace LKMC */
