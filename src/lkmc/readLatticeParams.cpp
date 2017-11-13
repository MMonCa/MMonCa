/*
 * readLatticeParam.cpp
 *
 *  Created on: Aug 16, 2012
 *
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

#include "io/ParameterManager.h"
#include "io/FileParameters.h"
#include "LatticeDiamondParam.h"
#include "LatticeBCC1Param.h"
#include "LatticeFCC1Param.h"
#include "LatticeDiamond.h"
#include "LatticeFCC1.h"
#include "LatticeBCC1.h"

using std::map;
using std::string;

namespace LKMC
{
	LatticeParam * readLatticeParams(Tcl_Interp *pTcl, const IO::ParameterManager *pPM, const IO::FileParameters *pPar, Kernel::M_TYPE mt)
	{
		if(pPM->getMaterialName(mt) == "Gas")
			return 0;
		std::string base = pPM->getMaterialName(mt) + "/Lattice/";
		std::string type = pPar->getString(base + "type");
		if(type == "diamond.1" || type == "diamond.2")
			return new LatticeDiamondParam(pTcl, pPM, pPar, mt);
		else if(type == "fcc.1")
			return new LatticeFCC1Param(pPM, pPar, mt);
		else if(type == "bcc.1")
			return new LatticeBCC1Param(pPM, pPar, mt);
		else
			ERRORMSG("Type '" << type << " not valid for lattice type");
		return 0;
	}

	Lattice * readLattice(Kernel::Domain *pD, const LatticeParam *p, Kernel::M_TYPE mt)
	{
		if(p == 0)
			return 0;
		switch(p->_type)
		{
		case LatticeParam::DIAMOND:
			return new LatticeDiamond(pD, p, mt);
		case LatticeParam::DIAMOND2:
			return new LatticeDiamond2(pD, p, mt);
		case LatticeParam::BCC:
			return new LatticeBCC1(pD, p, mt);
		case LatticeParam::FCC:
			return new LatticeFCC1(pD, p, mt);
		default:
			ERRORMSG("Lattice not implemented in ReadLattice");
			break;
		}
		return 0;
	}
}
