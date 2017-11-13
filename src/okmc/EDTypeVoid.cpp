/*
 * EDTypeVoid.cpp
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

#include "EDTypeVoid.h"
#include "Particle.h"
#include "io/ParameterManager.h"
#include "io/Diagnostic.h"
#include <algorithm>

using Kernel::Coordinates;
using std::map;


// disk, for instance a loop
namespace OKMC
{

std::vector<Kernel::Coordinates> EDTypeVoid::_template;

namespace temporal
{
	bool compare(const Kernel::Coordinates &i, const Kernel::Coordinates &j)
	{
		return i._x*i._x + i._y*i._y + i._z*i._z <  j._x*j._x + j._y*j._y + j._z*j._z;
	}
}

Coordinates EDTypeVoid::surface(Kernel::RNG &rng, const Kernel::Coordinates [2],
		const Kernel::Coordinates &center, const Kernel::ID &id, const Particle *pPart) const
{
	if(_template.size() == 0)
	{
		for(int i=-10; i<=10; ++i)
			for(int j=-10; j<=10; ++j)
				for(int k=-10; k<=10; ++k)
				{
					double a=1;
					_template.push_back(Coordinates(a*(i    ), a*(j    ), a*(k     )));
					_template.push_back(Coordinates(a*(i+.50), a*(j+.50), a*(k     )));
					_template.push_back(Coordinates(a*(i+.50), a*(j    ), a*(k+0.50)));
					_template.push_back(Coordinates(a*(i    ), a*(j+.50), a*(k+0.50)));
				}
		std::sort(_template.begin(), _template.end(), temporal::compare);
	}
	double templateDensity = 4.; //atoms/nm^3
	double scaling = std::pow(_densityNm / templateDensity, .3333333);
	unsigned size=0;
	for(map<Kernel::P_TYPE, unsigned>::const_iterator it=id._pt.begin(); it!=id._pt.end(); ++it)
		size += it->second;
	if(size > _template.size())
		ERRORMSG("Void is too big, maximum allowed size is " << _template.size());
	if(size == 0) size = 1;
	return _template[size-1]/scaling;
}

// sphere, for instance, a void
/*Coordinates EDTypeVoid::surface(Kernel::RNG &rng, const Kernel::Coordinates [2],
		const Kernel::Coordinates &center, const Kernel::ID &id, const Particle *pPart) const
{
	unsigned sizeImp = 0, sizeV = 0;
	map<P_TYPE, unsigned>::const_iterator it;
	for(it=id._pt.begin(); it!=id._pt.end(); ++it)
		if(it->first == V_TYPE)
			sizeV = it->second;
		else
			sizeImp += it->second;
	unsigned size;
	if(!pPart) //building from ID only
		size = std::max(sizeV, sizeImp);
	else
		size = (Domains::global()->PM()->getFamily(pPart->getPType()) == V_TYPE? sizeV : sizeImp);
	float r  = pow(size/_densityNm*3/(4*M_PI), .333333333); // 4/3 * pi * r^3 * density = size
	float theta = M_PI*rng.rand();
	float phi   = 2*M_PI*rng.rand();
	return Coordinates(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));

}*/

//tries to emit from the perifery (last incorporated particles)
const Particle * EDTypeVoid::emitFrom (Kernel::RNG &rng, const Kernel::Coordinates &center,
		const std::vector<Particle *> &parts) const
{
	if(parts.size() > 5)
	{
		unsigned idx = rng.rand()*5 + 1;
		return parts[parts.size() - idx];
	}
	return parts[rng.rand()*parts.size()];
}

}
