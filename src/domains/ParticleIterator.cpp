/*
 * ParticleIterator.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: ignacio.martin@imdea.org
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

#include "ParticleIterator.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "okmc/Particle.h"
#include "Global.h"

namespace Domains {

ParticleIterator::ParticleIterator() : _pPart(0)
{
}

ParticleIterator::ParticleIterator(const ParticleIterator &m)
{
	_pPart = m._pPart;
	_mei = m._mei;
}

ParticleIterator & ParticleIterator::operator=(const ParticleIterator &m)
{
	_pPart = m._pPart;
	_mei = m._mei;
	return *this;
}

ParticleIterator & ParticleIterator::operator++()
{
	if(_pPart)
		_pPart = _pPart->getNext();
	while(_pPart == 0)
	{
		++_mei;
		if(_mei == Domains::global()->endMEI())
		{
			_pPart =  0;
			break;
		}
		_pPart = _mei->getFirstPart();
	}
	return *this;
}

} /* namespace Domains */
