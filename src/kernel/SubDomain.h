/*
 * SubDomain.h
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

#ifndef SUBDOMAIN_H_
#define SUBDOMAIN_H_

#include "Selector.h"
#include "SelectorIterator.h"
#include "RNG.h"
#include "EventLog.h"
#include "ReactionLog.h"

namespace Kernel {
class Event;
class NullEvent;

class SubDomain {
public:
	typedef Selector<Event, double> Selector_f;

	SubDomain(Domain * p, unsigned idx, unsigned levels);
	~SubDomain();

	void insert(Event *, char, double);
	void update(Event *, char, double);
	void remove(Event *, char);
	void setkT(double);

	double getMaxRate() const;
	double nullEventRate(unsigned level);
	void   updateMaxRate(double maxrate, unsigned level);

	unsigned step(unsigned level);
	unsigned getIndex() const { return _idx; }
	unsigned getLevels() const { return _sel.size(); }

	Selector_f::iterator begin(char level) const { return _sel[level].begin(); }
	Selector_f::iterator end(char level) const { return _sel[level].end(); }

	Domain * _pDomain;
	RNG _rng;

	EventLog    _evLog;
	ReactionLog _reLog;
private:
	std::vector<Selector_f> _sel;
	std::vector<NullEvent *> _pNull;

	unsigned _idx;
	unsigned _levels;
};

} /* namespace Kernel */

#endif /* SUBDOMAIN_H_ */
