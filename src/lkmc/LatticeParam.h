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

#ifndef LKMCLATTICEPARAM_H
#define LKMCLATTICEPARAM_H

#include <tcl.h>
#include "kernel/Material.h"

namespace IO
{
	class FileParameters;
	class ParameterManager;
}
	
namespace LKMC {

	class LatticeParam
	{
	public:
		LatticeParam(const IO::ParameterManager *, const IO::FileParameters *, Kernel::M_TYPE);
		virtual ~LatticeParam() {}

		enum TYPE { DIAMOND, DIAMOND2, BCC, FCC };

		float _latticeParameter[2]; //[normal or alloy]
		TYPE _type             ;
		float _fIndex          [3];
		float _mIndex          [3];
		float _neighborDistance[3];
		float _growthDistortion;
		Kernel::M_TYPE _mt;

		double _orbitalRadius;
		bool   _mapToGrid;
	};

	LatticeParam * readLatticeParams(Tcl_Interp *, const IO::ParameterManager *, const IO::FileParameters *, Kernel::M_TYPE);
}

#endif
