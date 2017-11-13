/*
 * AlloyParam.h
 *
 *  Created on: Mar 11, 2013
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

#ifndef ALLOYPARAM_H_
#define ALLOYPARAM_H_
#include "kernel/Material.h"
#include "kernel/ParticleType.h"
#include "kernel/Event.h"
#include "kernel/Coordinates.h"
#include "io/ArrheniusAlloys.h"
#include "io/Polynomial.h"
#include <vector>


namespace IO { class FileParameters; class ParameterManager; }

namespace OKMC {

class AlloyParam
{
	static const unsigned SVAL = 3; // Number of neighbors considered for smoothing, 3 for Sides, Edges and, Corners
public:
	AlloyParam(const IO::ParameterManager *, const IO::FileParameters *);
	~AlloyParam() {}

	double getDerMixingEnergy(Kernel::M_TYPE, double, float); // Returns the x derivative of the energy value for the given composition
	double getMixingEnergy(Kernel::M_TYPE, double, float); // Returns the energy value for the given composition
        
    float  _corrFactor[Kernel::MAX_MATERIALS][2]; //Self-diffusion correlation factors (fi, fv)
	float  _smoothing[Kernel::MAX_MATERIALS][SVAL];  //Smoothing factors for [0] Side boxes, [1] Edge boxes, [2] Corner boxes
	IO::ArrheniusAlloys _alpha[Kernel::MAX_MATERIALS][2]; // Self diffusion coefficient for the AB alloy (alpha = D_B / D_A)
	
	float _theta[Kernel::MAX_MATERIALS];
    bool _isAlloy[Kernel::MAX_MATERIALS]; // Returns true if the simulated material is an alloy, false otherwise

    IO::Polynomial _mixingEnthalpy[Kernel::MAX_MATERIALS];

#ifdef NUMODEL
    double nu(double x, double kT);
#endif
};
}

#endif /* ALLOYPARAM_H_ */
