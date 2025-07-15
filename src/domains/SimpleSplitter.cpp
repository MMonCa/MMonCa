/*
 * SimpleSplitter.cpp
 *
 *  Created on: Aug 6, 2012
 *      Author: ignacio.martin@imdea.org
 *
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

#include "SimpleSplitter.h"
#include "Global.h"
#include "io/FileParameters.h"
#include <limits>

namespace Domains {

SimpleSplitter::SimpleSplitter(
		const Kernel::Coordinates &mm, const Kernel::Coordinates &MM,  //global simulation size
		const IO::GetMaterial *gm  //pointer to getMaterial function
		) : Splitter(mm, MM, gm)
{
	_nDomains = global()->getFileParameters()->getInt("MC/General/domains");
	unsigned subDomains = global()->getFileParameters()->getInt("MC/General/subdomains");
	_zSlide = (_maxCell._z - _minCell._z)/_nDomains;
	MEDMSG(this << " Using " << _nDomains << " domains with a z of " << _zSlide);

	for(unsigned dom=0; dom < _nDomains; ++dom)
	{
		_subDomains.push_back(std::vector<box>());
		Kernel::Coordinates m(mm), M(MM);
		splitDomain(dom, m, M);
		switch(subDomains)
		{
		case 1:
			_subDomains.back().push_back(box(m, M));
			break;
		case 2:
			_subDomains.back().push_back(box(m, M));
			_subDomains.back().push_back(box(m, M));
			_subDomains.back()[0]._max._z = (m._z + M._z)/2;
			_subDomains.back()[1]._min._z = (m._z + M._z)/2;
			break;
		case 4:
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back()[0]._max._z = (m._z + M._z)/2;
			_subDomains.back()[0]._max._y = (m._y + M._y)/2;
			_subDomains.back()[1]._max._z = (m._z + M._z)/2;
			_subDomains.back()[1]._min._y = (m._y + M._y)/2;
			_subDomains.back()[2]._min._z = (m._z + M._z)/2;
			_subDomains.back()[2]._max._y = (m._y + M._y)/2;
			_subDomains.back()[3]._min._z = (m._z + M._z)/2;
			_subDomains.back()[3]._min._y = (m._y + M._y)/2;
			break;
		case 6:
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back()[0]._max._y = m._y + (M._y - m._y)/3.;
			_subDomains.back()[1]._min._y = _subDomains.back()[0]._max._y;
			_subDomains.back()[1]._max._y = m._y + 2*(M._y - m._y)/3.;
			_subDomains.back()[2]._min._y = _subDomains.back()[1]._max._y;
			_subDomains.back()[3]._max._y = m._y + (M._y - m._y)/3.;
			_subDomains.back()[4]._min._y = _subDomains.back()[3]._max._y;
			_subDomains.back()[4]._max._y = m._y + 2*(M._y - m._y)/3.;
			_subDomains.back()[5]._min._y = _subDomains.back()[4]._max._y;
			_subDomains.back()[0]._max._z = (m._z + M._z)/2;
			_subDomains.back()[1]._max._z = (m._z + M._z)/2;
			_subDomains.back()[2]._max._z = (m._z + M._z)/2;
			_subDomains.back()[3]._min._z = (m._z + M._z)/2;
			_subDomains.back()[4]._min._z = (m._z + M._z)/2;
			_subDomains.back()[5]._min._z = (m._z + M._z)/2;
			break;
		case 8:
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));
			_subDomains.back().push_back(box(m,M));

			_subDomains.back()[0]._max._y = m._y + (M._y - m._y)/4.;
			_subDomains.back()[1]._min._y = _subDomains.back()[0]._max._y;
			_subDomains.back()[1]._max._y = m._y + 2*(M._y - m._y)/4.;
			_subDomains.back()[2]._min._y = _subDomains.back()[1]._max._y;
			_subDomains.back()[2]._max._y = m._y + 3*(M._y - m._y)/4.;
			_subDomains.back()[3]._min._y = _subDomains.back()[2]._max._y;

			_subDomains.back()[4]._max._y = _subDomains.back()[0]._max._y;
			_subDomains.back()[5]._min._y = _subDomains.back()[1]._min._y;
			_subDomains.back()[5]._max._y = _subDomains.back()[1]._max._y;
			_subDomains.back()[6]._min._y = _subDomains.back()[2]._min._y;
			_subDomains.back()[6]._max._y = _subDomains.back()[2]._max._y;
			_subDomains.back()[7]._min._y = _subDomains.back()[3]._min._y;

			_subDomains.back()[0]._max._z = (m._z + M._z)/2;
			_subDomains.back()[1]._max._z = (m._z + M._z)/2;
			_subDomains.back()[2]._max._z = (m._z + M._z)/2;
			_subDomains.back()[3]._max._z = (m._z + M._z)/2;

			_subDomains.back()[4]._min._z = (m._z + M._z)/2;
			_subDomains.back()[5]._min._z = (m._z + M._z)/2;
			_subDomains.back()[6]._min._z = (m._z + M._z)/2;
			_subDomains.back()[7]._min._z = (m._z + M._z)/2;
			break;
		default:
			ERRORMSG("Sorry, but I cannot accept " << subDomains << " subdomains.");
			break;
		}
	}
}


void SimpleSplitter::splitDomain(unsigned nDomain, Kernel::Coordinates &m, Kernel::Coordinates &M) const //given a domain number, returns the block it occupies
{
	m = _minCell;
	M = _maxCell;
	m._z = _minCell._z + nDomain*_zSlide;
	M._z = _minCell._z + (nDomain+1)*_zSlide;
}

std::vector<float> SimpleSplitter::getLinesZpart(float const aZmin, float const aZmax, std::vector<float> const * const aLinesZ) {
	uint32_t indexMin = getClosestIndex(aZmin, aLinesZ);	
	uint32_t indexMax = getClosestIndex(aZmax, aLinesZ);
	std::vector<float> result;
	for(uint32_t i = indexMin; i <= indexMax; ++i) {
		result.push_back(aLinesZ->at(i));
	}
	return result;
}

unsigned SimpleSplitter::getDomain(const Kernel::Coordinates &c) const
{
	double f = (double(c._z) - double(_minCell._z)) / double(_zSlide);
	return f;
}

unsigned SimpleSplitter::getSubDomain(const Kernel::Coordinates &m, const Kernel::Coordinates &M) const
{
	Kernel::Coordinates center = m;
	for(unsigned dom=0; dom<_nDomains; ++dom)
		for(unsigned i=0; i<_subDomains[dom].size(); ++i)
			if(center.isInto(_subDomains[dom][i]._min, _subDomains[dom][i]._max))
				return i;
	ERRORMSG("Cannot find subdomain");
	return 0;
}

unsigned SimpleSplitter::getLevel(const Kernel::Coordinates &m, const Kernel::Coordinates &M) const
{
	if(_levels == 1)
		return 0;
	if(_levels != 9 && _levels != 3)
		ERRORMSG("levels can only be 1, 3 or 9");
	if(_levels == 3)
		switch(_subDomains.size())
		{
		case 4: case 6: case 8:
			ERRORMSG("Sorry, but 9 levels are compulsory with 4 " << _subDomains.size() << " subdomains.");
			break;
		};


	Kernel::Coordinates minCell, maxCell, center;
	unsigned sub = getSubDomain(m, M);
	center = m;
	unsigned dom = getDomain(m);
	getSubDomain(dom, sub, minCell, maxCell);

	double slidey = (maxCell._y - minCell._y)/3.;
	double slidez = (maxCell._z - minCell._z)/3.;
	unsigned y = (center._y - minCell._y) / slidey;
	unsigned z = (center._z - minCell._z) / slidez;

	if(_levels == 3)
		return z;
	return y*3+z;
}

uint32_t SimpleSplitter::getClosestIndex(float const aWhere, std::vector<float> const * const aLines) {
	uint32_t result = 0u;
	float min = std::numeric_limits<float>::max();
	for(uint32_t i = 0u; i < aLines->size(); ++i) { 
		auto const diff = std::abs(aWhere - aLines->at(i));
		if(diff < min) {
			min = diff;
			result = i;
		}
	}
	return result;
}

} /* namespace Domains */
