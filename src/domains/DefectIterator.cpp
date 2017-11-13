/*
 * DefectIterator.cpp
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

#include "DefectIterator.h"
#include "kernel/Selector.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/RateManager.h"
#include "domains/Splitter.h"
#include "okmc/Defect.h"

namespace Domains {


DefectIterator::DefectIterator(const Kernel::SelectorIterator<Kernel::Event, double> &si) :
		_pSi(new Kernel::SelectorIterator<Kernel::Event, double>(si))
{
	_nDomain = 0; _nSubDomain = 0; _nLevel = 0;
}

DefectIterator::DefectIterator(const DefectIterator &di) : _pSi(new Kernel::SelectorIterator<Kernel::Event, double>(*di._pSi))
{
	_nDomain    = di._nDomain;
	_nSubDomain = di._nSubDomain;
	_nLevel     = di._nLevel;
}

DefectIterator & DefectIterator::operator=(const DefectIterator &di)
{
	_pSi = di._pSi;
	_nDomain    = di._nDomain;
	_nSubDomain = di._nSubDomain;
	_nLevel     = di._nLevel;
	return *this;
}

DefectIterator::~DefectIterator()
{
}

DefectIterator & DefectIterator::operator++()
{
	do
	{
		if(*_pSi != global()->getDomain(_nDomain)->_pRM->getSubDomain(_nSubDomain)->end(_nLevel))
			++(*_pSi);
		while( *_pSi == global()->getDomain(_nDomain)->_pRM->getSubDomain(_nSubDomain)->end(_nLevel)) //try next level
		{
			_nLevel++;
			if(_nLevel == global()->getSplitter()->getLevels()) //try next subdomain
			{
				_nSubDomain++;
				_nLevel = 0;
				if(_nSubDomain == global()->getSplitter()->getSubDomains(_nDomain)) //try next domain
				{
					_nDomain++;
					_nSubDomain = 0;
					if(_nDomain == global()->getDomains())
					{
						_pSi.reset();
						return *this;
					}
				}
			}
			*_pSi = global()->getDomain(_nDomain)->_pRM->getSubDomain(_nSubDomain)->begin(_nLevel);
		}
	} while(dynamic_cast<OKMC::Defect *>(**_pSi) == 0);

	return *this;
}

bool DefectIterator::operator!=(const DefectIterator &m) const
{
	if(m._pSi == 0 && _pSi == 0)
		return false;
	if(m._pSi == 0 || _pSi == 0)
		return true;
	return *m._pSi != *_pSi;
}

const OKMC::Defect * DefectIterator::operator*()  const
{
	return dynamic_cast<OKMC::Defect *>(**_pSi);
}

const OKMC::Defect * DefectIterator::operator->() const
{
	return dynamic_cast<OKMC::Defect *>(**_pSi);
}



} /* namespace Domains */
