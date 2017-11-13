/*
 * MeshParam.cpp
 *
 *  Created on: Feb 10, 2011
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

#include "MeshParam.h"
#include "io/FileParameters.h"
#include "io/ParameterManager.h"

namespace Kernel {

MeshParam::MeshParam(const IO::ParameterManager *pPM, const IO::FileParameters *pPar)
{
	_periodicX  = pPar->getBool ("MC/Mesh/periodic.x");
	_periodicY  = pPar->getBool ("MC/Mesh/periodic.y");
	_periodicZ  = pPar->getBool ("MC/Mesh/periodic.z");

	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
		_lambda[mt] = pPar->getFloat(pPM->getMaterialName(mt) + "/Models/lambda");

	_bLHF       = pPar->getBool ("MC/Mesh/speed.up");

	_spacing[0] = pPar->getFloat("MC/Mesh/spacing.x");
	_spacing[1] = pPar->getFloat("MC/Mesh/spacing.y");
	_spacing[2] = pPar->getFloat("MC/Mesh/spacing.z");

	_poissonX   = pPar->getString("MC/Mesh/poisson.bc.x");
	_poissonY   = pPar->getString("MC/Mesh/poisson.bc.y");
	_poissonZ   = pPar->getString("MC/Mesh/poisson.bc.z");

	// value of the potential at x, y or z edge
	// TODO: instead of defining the potential for an edge, it should be defined for a node (use a Tcl procedure)
	_potentialX = pPar->getFloat("MC/Mesh/potential.x");
	_potentialY = pPar->getFloat("MC/Mesh/potential.y");
	_potentialZ = pPar->getFloat("MC/Mesh/potential.z");
}

MeshParam::~MeshParam() { }

}
