/*
 * NullEvent.h
 *
 *  Created on: Oct 7, 2013
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

#ifndef NULLEVENT_H_
#define NULLEVENT_H_

#include "Event.h"

namespace Kernel {

class NullEvent: public Kernel::Event {
public:
	NullEvent(Kernel::Domain *pD, char l) : Event(pD), _rate(0), _level(l) {}
	virtual ~NullEvent() {}

	void                perform (Kernel::SubDomain *, unsigned );
	void                setIndex(unsigned, int idx) {_idx = idx; }
	unsigned            getNEvents() const { return 1; }
	int                 getIndex(unsigned ) const { return _idx; }
	float               getRate(unsigned , float ) const { return _rate; }
	char                getSubDomainLevel() const { return _level; }
	E_TYPE              getEType() const { return Event::EMPTY; }

	double _rate;

private:
	int _idx;
	char _level;
};

} /* namespace Kernel */

#endif /* NULLEVENT_H_ */
