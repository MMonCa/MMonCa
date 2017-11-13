/*
 * EDType.cpp
 *
 *  Created on: Jul 2, 2013
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

#include "EDType.h"
#include "io/Diagnostic.h"
using namespace OKMC;
using Kernel::MAX_PARTICLES;
using Kernel::MAX_IMPURITIES;
using Kernel::ID;
using Kernel::P_TYPE;
using std::map;

EDType::CLType::CLType() {
	for (int i = 0; i < MAX_PARTICLES; ++i) {
		_pref[i] = 0;
		_interactMP[i] = false;
		_interactMP_sink[i] = false;
		_interactMP_radius[i] = 0;
	}
}

void EDType::checkAxes()
{
	if(_axes[0] * _axes[1] != 0)
		ERRORMSG("In defect " << _name << " axis.0 not perpendicular to axis.1");
	_axes[2] = _axes[0].product(_axes[1]);
	for(int i=0; i<3; ++i)
		_axes[i] /= _axes[i].abs();
}

