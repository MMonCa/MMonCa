/*
 * SubDomain.cpp
 *
 *  Created on: Oct 7, 2013
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

#include "SubDomain.h"
#include "NullEvent.h"
#include "Domain.h"
#include "RateManager.h"
#include "okmc/ClusterParam.h"
#include <tcl.h>

namespace Kernel {

SubDomain::SubDomain(Domain *p, unsigned idx, unsigned levels) : _pDomain(p),
		_rng(_pDomain->_rng_dom.rand()*5000), _idx(idx), _levels(levels)
{
	for(unsigned level=0; level<_levels; level++)
	{
		_sel.push_back(Selector_f());
		_pNull.push_back(new NullEvent(_pDomain, level));
		unsigned idx = _sel[level].insert(_pNull[level], 0., 0);
		_pNull.back()->setIndex(0, idx);
	}
}

SubDomain::~SubDomain()
{
	for(unsigned level=0; level<_levels; level++)
	{
		for(Selector_f::iterator it = _sel[level].begin(); it != _sel[level].end(); ++it)
			delete *it;
		_pNull[level] = 0;
	}
}

unsigned SubDomain::step(unsigned level)
{
	int tag;
	Event * pEvent;
	_sel[level].select(_rng.rand(), pEvent, tag);
	pEvent->perform(this, tag);
	return (pEvent == _pNull[level]? 0:1);
}

double SubDomain::getMaxRate() const
{
	double maxRate = 0;
	for(unsigned i=0; i<_levels; ++i)
	{
		double realRate = _sel[i].getTotalRate() - _pNull[i]->_rate;
		if(realRate > maxRate)
			maxRate = realRate;
	}
	return maxRate;
}


double SubDomain::nullEventRate(unsigned level)
{
	_pNull[level]->_rate = 0;
	_sel[level].update(_pNull[level]->getIndex(0), _pNull[level]->_rate);
	return _sel[level].getTotalRate();
}


void SubDomain::updateMaxRate(double maxrate, unsigned level)
{
	double realRate = _sel[level].getTotalRate() - _pNull[level]->_rate;
	_pNull[level]->_rate = maxrate - realRate;
	if(_pNull[level]->_rate < 0)
	{
		ERRORMSG(level << " Setting negative null rate total=" << "( " << maxrate << " - " << realRate << " ) = " << _pNull[level]->_rate << " to zero");
		_pNull[level]->_rate = 0;
	}
	_sel[level].update(_pNull[level]->getIndex(0), _pNull[level]->_rate);
}

void SubDomain::setkT(double kT)
{
	for(unsigned level=0; level <_levels; ++level)
		for(Selector_f::iterator it = begin(level); it != end(level); ++it)
			update(*it, level, kT);
}

void SubDomain::update(Event *mA, char level, double kT)
{
	for(unsigned ev=0; ev < mA->getNEvents(); ++ev)
		_sel[level].update(mA->getIndex(ev), mA->getRate(ev, kT));
}

void SubDomain::insert(Event *pEv, char level, double kT)
{
	for(unsigned ev=0; ev < pEv->getNEvents(); ++ev)
	{
		unsigned idx = _sel[level].insert(pEv, pEv->getRate(ev, kT), ev);
		pEv->setIndex(ev, idx);
	}
}

void SubDomain::remove(Event *pEv, char level)
{
	for(int ev=pEv->getNEvents()-1; ev >=0; --ev)
		_sel[level].remove(pEv->getIndex(ev));
}

} /* namespace Kernel */
