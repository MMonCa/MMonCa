/* Author: ignacio.martin@imdea.org
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

#include "RateManager.h"
#include "Domain.h"

#include "io/FileParameters.h"
#include "domains/Splitter.h"
#include "SubDomain.h"
#include "MeshElement.h"
#include "electrostatics/Poisson.h"
#include "mechanics/MechInterface.h"
#include <iostream>
#include <limits>
#include <iomanip>

using LKMC::LatticeAtom;
using std::string;
using std::vector;

namespace Kernel {

RateManager::RateManager(Domain *p) :
		_pDomain(p),
		_kB(8.6174e-5)//in eV/K
{
	_time        = _timeLastEvent = _timeNextEvent = 0;
	_lastMaxRate = 0;
	_nEvents     = 0;
	_kelvin      = 23 + 273.15;
	_kT          = _kB * _kelvin;
	_depthLA     = std::numeric_limits<float>::max();

	unsigned nSubDomains = Domains::global()->getSplitter()->getSubDomains(_pDomain->_domain_number);
	_nLevels     = Domains::global()->getSplitter()->getLevels();
	_bAveraged   = Domains::global()->getFileParameters()->getBool("MC/General/average.rates");
	for(unsigned i=0; i<nSubDomains; ++i)
		_subDomains.push_back(new SubDomain(_pDomain, i, _nLevels));
}

RateManager::~RateManager()
{
	while(_subDomains.size())
	{
		delete _subDomains.back();
		_subDomains.pop_back();
	}
}

void RateManager::setTempK(float K)
{
	if(K == _kelvin)
		return;
	_kelvin = K;
	_kT     = _kB * K;

	if(Domains::global()->IsPoisson())
		_pDomain->_pPoisson->setT(_kelvin);

	for(vector<SubDomain *>::iterator i = _subDomains.begin(); i != _subDomains.end(); ++i)
		(*i)->setkT(_kT);
}

void RateManager::update(Event *pEv, MeshElement *pEle)
{
	char idx = pEle->getSubDomainIdx();
	_subDomains[idx]->update(pEv, pEle->getLevel(), _kT);
}

void RateManager::insert(Event *pEv, MeshElement *pEle)
{
	char idx = pEle->getSubDomainIdx();
	_subDomains[idx]->insert(pEv, pEle->getLevel(), _kT);
}

void RateManager::remove(Event *pEv, MeshElement *pEle)
{
	char idx = pEle->getSubDomainIdx();
	_subDomains[idx]->remove(pEv, pEle->getLevel());
}

//annealing time OR LKMC depth in nanometers
void RateManager::anneal(double time, bool bDepth, float depth, long unsigned events)
{
	double endTime = _time+time;
	while(_timeNextEvent < endTime)
	{
		double maxRate = 0;
		for(vector<SubDomain *>::iterator it = _subDomains.begin(); it != _subDomains.end(); ++it)
		{
			double gmr = (*it)->getMaxRate();
			if(maxRate < gmr)
				maxRate = gmr;
		}
		double rand, deltaTime = 0;

		if(maxRate != 0)
		{
			while ( (rand = _pDomain->_rng_dom.rand()) == 0);
			deltaTime += (_bAveraged? 1. : -std::log(rand)) / (_nLevels*maxRate);
			_lastMaxRate   = maxRate;
			_timeLastEvent = _timeNextEvent;
			_timeNextEvent = std::min(endTime, _timeNextEvent + deltaTime);
			_time          = _timeNextEvent;
		}
		else
		{
			WARNINGMSG("No more events to anneal... over.");
			LOWMSG("No more events to anneal... over.");
			_time = endTime;
			break;
		}
		unsigned level = (_nLevels == 1? 0 : (_pDomain->_rng_dom.rand()*_nLevels));
		unsigned valid = 0;
#ifdef _RELEASE_PARALLEL
		#pragma omp parallel num_threads(_subDomains.size()) reduction(+:valid)
		{
			#pragma omp for schedule(static,1)
#endif
			for(unsigned nd = 0; nd < _subDomains.size(); ++nd)
			{
				_subDomains[nd]->updateMaxRate(maxRate, level);
				valid += _subDomains[nd]->step(level);
			}
#ifdef _RELEASE_PARALLEL
		}
#endif
		_nEvents += valid;

		if(Domains::global()->mechModel() != "None" && _pDomain->_um_mechani(_time, _nEvents, _depthLA))
		{
			_pDomain->_pMech->import();
			for(vector<SubDomain *>::iterator i = _subDomains.begin(); i != _subDomains.end(); ++i)
				(*i)->setkT(_kT);
		}
		if(Domains::global()->IsPoisson() && _pDomain->_um_poisson(_time, _nEvents, _depthLA))
			_pDomain->_pPoisson->compute();

		if( (bDepth && _depthLA < depth) || (events > 0 && _nEvents >= events)
				|| (Domains::global()->getDomains() == 1 && _subDomains.size() == 1 && Domains::global()->getNextEvents() < _nEvents) )
			break;
	}
}

void RateManager::setDepthLA(float depth)
{
	if(depth < _depthLA)
		_depthLA = depth;
}

void RateManager::restart(std::ostream &os) const
{
	os <<  _nEvents << " " << _depthLA << " " <<
			_time << " " << _timeLastEvent << " " <<
			_timeNextEvent << " " << _lastMaxRate;
}

void RateManager::restart(std::istream &is)
{
	is >>  _nEvents >> _depthLA >> _time >> _timeLastEvent >> _timeNextEvent >> _lastMaxRate;
}

}
