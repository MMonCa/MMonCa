/*
 * EDTypeDisk.cpp
 *
 *  Created on: Aug 14, 2013
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

#include "EDTypeDisk.h"
#include "Particle.h"
#include "io/ParameterManager.h"
#include "domains/Global.h"

using Kernel::Coordinates;
using std::map;

// disk, for instance a loop
namespace OKMC
{

/* For I particles, it builds a disk
 * For all other particles, they are left where they were found.
 * If no particles are specified, it builds a disk.
 */
Coordinates EDTypeDisk::surface(Kernel::RNG & rng, const Kernel::Coordinates c[2], const Kernel::Coordinates &center,
		const Kernel::ID &id, const Particle *pPart) const
{
	const float l = 0.3;
	if(_intState.size() == 0) //init intSta. For every layer it has how many defects can manage.
	{
		_intState.push_back(0);
		for(unsigned i=2; i<2000; ++i)
		{
			unsigned int n = ceil(sqrt(i/(M_PI*_densityNm))/l); //"radius quantization"
			if(n > _intState.size())
				_intState.push_back(i-1);
		}
	}

	map<Kernel::P_POS, unsigned>::const_iterator it = id._pos.find(Kernel::POS_I);
	unsigned size = (it == id._pos.end()? 1 : it->second);

	if(pPart == 0 || Domains::global()->PM()->isIorV(id._mt, pPart->getPType()) == 0)
	{
		int n = ceil(sqrt(size/(M_PI*_densityNm))/l); //"radius quantization"
		float radius = n*l;
		int howManyHere = _intState[n] - _intState[n-1];
		float alpha = 2*M_PI*float(size-_intState[n-1])/float(howManyHere);
		float Ca = cos(alpha);
		float Sa = sin(alpha);
		return c[0]*radius*Ca + c[1]*radius*Sa;
	}
	else
		return pPart->getCoordinates() - center;
}

const Particle * EDTypeDisk::emitFrom (Kernel::RNG & rng, const Kernel::Coordinates &center,
				const std::vector<Particle *> &parts) const
{
	return parts[rng.rand()*parts.size()];
}
}

