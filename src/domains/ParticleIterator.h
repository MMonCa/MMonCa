/*
 * ParticleIterator.h
 *
 *  Created on: Jun 20, 2012
 *      Author: ignacio.martin@imdea.org
 *
 *      Iterator to obtain particles
 *      It is "pseudo", because it does not contain a pointer to a value
 *      but it is, somehow, the value itself.
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

#ifndef PARTICLEITERATOR_H_
#define PARTICLEITERATOR_H_

#include "MeshElementIterator.h"

namespace OKMC   { class Particle; }

namespace Domains {

class ParticleIterator {
public:
	ParticleIterator();
	ParticleIterator(const ParticleIterator &m);
	ParticleIterator & operator=(const ParticleIterator &m);

	ParticleIterator & operator++();
	bool               operator!=(const ParticleIterator &m) const { return m._pPart != _pPart; }
	bool               operator==(const ParticleIterator &m) const { return !operator!=(m); }

	const OKMC::Particle * operator*()  const { return _pPart; }
	const OKMC::Particle * operator->() const { return _pPart; }
	OKMC::Particle *       modify()           { return _pPart; }

private:
	MeshElementIterator _mei;
	OKMC::Particle *_pPart;

	friend class Domains::Global;
};

} /* namespace Domains */

#endif /* PARTICLEITERATOR_H_ */
