/*
 * MeshParam.h
 *
 *  Created on: Feb 10, 2011
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

#ifndef MESHPARAM_H_
#define MESHPARAM_H_

#include <string>
#include "Material.h"

namespace IO
{
	class FileParameters;
	class ParameterManager;
}

namespace Kernel {

class MeshParam {
public:
	MeshParam(const IO::ParameterManager *, const IO::FileParameters *);
	virtual ~MeshParam();

	bool        _periodicX, _periodicY, _periodicZ;
	bool        _bLHF; //long hops on/off
	std::string _poissonX, _poissonY, _poissonZ;
	double      _potentialX, _potentialY, _potentialZ;

	float _lambda[Kernel::MAX_MATERIALS];
	float _spacing[3];
};

}

#endif /* MESHPARAM_H_ */
