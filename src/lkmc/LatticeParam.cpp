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

#include "LatticeParam.h"
#include "io/ParameterManager.h"
#include "io/FileParameters.h"

using std::string;
using std::map;

namespace LKMC {

LatticeParam::LatticeParam(const IO::ParameterManager *pPM, const IO::FileParameters *pPar, Kernel::M_TYPE mt) : _mt(mt)
{
	string base =  pPM->getMaterialName(mt) + "/Lattice/";

	_latticeParameter[0] = pPar->getFloat (base + "parameter");

	if(pPM->getMaterial(mt)._alloy != Kernel::UNDEFINED_TYPE)
		_latticeParameter[1] = pPar->getFloat (base + "parameter.alloy");
	else
		_latticeParameter[1] = _latticeParameter[0];
		
	map<string, float> floatMap = pPar->getFloatMap(base + "wafer.orientation");
	if( floatMap.find("i") == floatMap.end() ||
			floatMap.find("j") == floatMap.end() ||
			floatMap.find("k") == floatMap.end())
		ERRORMSG("indexes i, j, k not found in wafer.orientation");
	_mIndex[0] = floatMap["i"];
	_mIndex[1] = floatMap["j"];
	_mIndex[2] = floatMap["k"];
	
	floatMap = pPar->getFloatMap(base + "flat.orientation");
	if( floatMap.find("i") == floatMap.end() ||
			floatMap.find("j") == floatMap.end() ||
			floatMap.find("k") == floatMap.end())
		ERRORMSG("indexes i, j, k not found in flat.orientation");
	_fIndex[0] = floatMap["i"];
	_fIndex[1] = floatMap["j"];
	_fIndex[2] = floatMap["k"];
	
	floatMap = pPar->getFloatMap(base + "neighbor.distance");
	if(floatMap.find("first") == floatMap.end() ||
	  floatMap.find("second") == floatMap.end() ||
	  floatMap.find("third") == floatMap.end())
	  ERRORMSG("first and second should be indexes of " << base + "neighbor.distance");
	_neighborDistance[0] = floatMap["first"];
	_neighborDistance[1] = floatMap["second"];
	_neighborDistance[2] = floatMap["third"];

	if (pPar->specified(base + "map.to.grid") && pPar->specified(base + "orbital.radius")) {
		_mapToGrid     = pPar->getBool(base + "map.to.grid");
		_orbitalRadius = pPar->getFloat(base + "orbital.radius");
	}
	else {
		_mapToGrid     = false;
		_orbitalRadius = 0.;
	}

	if (pPar->specified(base+"growth.distortion"))
		_growthDistortion = pPar->getFloat(base + "growth.distortion");
	else
		_growthDistortion = 0;
}

}
