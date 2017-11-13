/*
 * EDTypes.cpp
 *
 *  Created on: Jun 29, 2011
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

#include "EDTypeIrregular.h"
#include "Particle.h"

using Kernel::Coordinates;

namespace OKMC {

/*
 *  Irregular, for instance a small cluster
 *  If no particles specified, it creates an irregular sphere.
 */

Coordinates EDTypeIrregular::surface(Kernel::RNG &rng, const Kernel::Coordinates [2],
		const Kernel::Coordinates &center, const Kernel::ID &id, const Particle *pPart) const
{
	if(pPart)
		return pPart->getCoordinates() - center;
	unsigned size = 0;
	for(std::map<Kernel::P_TYPE, unsigned>::const_iterator it = id._pt.begin(); it != id._pt.end(); ++it)
		size += it->second;
	double radius = pow(3*size/(4*M_PI*_densityNm), 1./3.);
	double theta = rng.rand() * 2*M_PI;
	double phi   = rng.rand() * 2*M_PI;
	return Coordinates(radius*cos(theta)*sin(phi), radius*sin(theta), radius*cos(theta)*cos(phi));
}

const Particle * EDTypeIrregular::emitFrom (Kernel::RNG & rng, const Kernel::Coordinates &center,
		const std::vector<Particle *> &parts) const
{
	return parts.back();
}

}
