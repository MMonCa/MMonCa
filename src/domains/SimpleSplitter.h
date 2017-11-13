/*
 * SimpleSplitter.h
 *
 *  Created on: Aug 6, 2012
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

#ifndef SIMPLESPLITTER_H_
#define SIMPLESPLITTER_H_

#include "Splitter.h"
#include <vector>

namespace Domains {

class SimpleSplitter: public Splitter {
public:
	SimpleSplitter(const Kernel::Coordinates &m, const Kernel::Coordinates &M,  //global simulation size
			const IO::GetMaterial *  //pointer to getMaterial function
			);
	virtual ~SimpleSplitter() {}

	virtual unsigned getDomains() const { return _nDomains; }
	virtual unsigned getSubDomains(unsigned dom) const { return _subDomains[dom].size(); }
	virtual void splitDomain(unsigned nDomain, Kernel::Coordinates &m, Kernel::Coordinates &M) const;
	virtual unsigned getDomain(const Kernel::Coordinates &) const;
	virtual unsigned getSubDomain(const Kernel::Coordinates &m, const Kernel::Coordinates &M) const;
	virtual unsigned getLevel(const Kernel::Coordinates &m, const Kernel::Coordinates &M) const;

private:
	struct box {
		box(const Kernel::Coordinates &m, const Kernel::Coordinates &M) : _min(m), _max(M) {}
		Kernel::Coordinates _min;
		Kernel::Coordinates _max;
	};
	void getSubDomain(unsigned dom, unsigned i, Kernel::Coordinates &m, Kernel::Coordinates &M) const
		{ m=_subDomains[dom][i]._min; M=_subDomains[dom][i]._max; }
	std::vector<std::vector<box> > _subDomains; //domains, subdomains.
	unsigned _nDomains;
	float _zSlide;
};

} /* namespace Domains */
#endif /* SIMPLESPLITTER_H_ */
