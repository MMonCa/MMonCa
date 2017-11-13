/*
 * Splitter.h
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

#ifndef SPLITTER_H_
#define SPLITTER_H_

#include "kernel/Coordinates.h"

namespace IO { class GetMaterial; }

namespace Domains {

class Splitter {
public:
	Splitter(const Kernel::Coordinates &m, const Kernel::Coordinates &M, const IO::GetMaterial *);
	virtual ~Splitter() {}

	virtual unsigned getDomains()    const = 0;  //returns how many domains are here
	virtual unsigned getSubDomains(unsigned nDom) const = 0;
	        unsigned getLevels()     const { return _levels; }
	virtual void splitDomain(unsigned nDomain, Kernel::Coordinates &m, Kernel::Coordinates &M) const = 0; //given a domain number, returns the block it occupies
	virtual unsigned getDomain(const Kernel::Coordinates &) const = 0; //given a global coordinate, returns the domain.
	virtual unsigned getSubDomain(const Kernel::Coordinates &m, const Kernel::Coordinates &M) const = 0;
	virtual unsigned getLevel(const Kernel::Coordinates &m, const Kernel::Coordinates &M) const = 0;

protected:
	Kernel::Coordinates _minCell;
	Kernel::Coordinates _maxCell;  //global simulation size
	const IO::GetMaterial *_getMaterial; //pointer to getMaterial function
	char     _levels;
};

} /* namespace Kernel */
#endif /* SPLITTER_H_ */
