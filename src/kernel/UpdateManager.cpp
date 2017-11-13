/*
 * UpdateManager.cpp
 *
 *  Created on: Nov 7, 2013
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

#include "UpdateManager.h"
#include "domains/Global.h"
#include "io/FileParameters.h"
#include <limits>

using std::string;

namespace Kernel {

UpdateManager::UpdateManager(const string &tag) : _tag(tag)
{
	_bPrint = true;
	_decadesT = Domains::global()->getFileParameters()->getInt(tag+".time.decade");
	_deltaT   = Domains::global()->getFileParameters()->getFloat(tag+".time.delta");
	_howMany  = Domains::global()->getFileParameters()->getFloat(tag+".events");
	_minT     = Domains::global()->getFileParameters()->getFloat(tag+".time.min");
	_oldEvents = 0;
	if(_minT == 0)
		_minT = 1e-20;
	_nextTime = _minT;
	_nextEvents = _howMany;
}

bool UpdateManager::operator()(double time, unsigned long int events, double depth) const
{
	bool ret = false;
	while((_deltaT != 0 || _decadesT != 0) && time >= _nextTime)
	{
		_nextTime = computeNextTime();
		ret = true;
	}
	while(_howMany != 0 && events >= _nextEvents)
	{
		_nextEvents = computeNextEvent();
		ret = true;
	}
	if(ret && _bPrint)
		LOWMSG("Updating " << _tag << " Next: " << _nextTime << " - " << _nextEvents);
	if(_oldEvents != 0 && _oldEvents == events)
		return false;
	_oldEvents = events;
	return ret;
}

double UpdateManager::computeNextTime() const
{
	double nextTimeDelta = _deltaT ? (unsigned(_nextTime / _deltaT) + 1)*_deltaT : std::numeric_limits<double>::max();
	double nextTimeDecad = std::numeric_limits<double>::max();
	if(_decadesT)
	{
		int p10 = int(log10(_nextTime));
		if(p10 < 0)
			--p10;
		do
		{
			double dT = 10./_decadesT;
			dT *= pow(10., p10);
			nextTimeDecad = (unsigned(_nextTime/dT) +1)*dT;
			p10++;
		} while(nextTimeDecad == _nextTime);
	}
	return (nextTimeDelta < nextTimeDecad? nextTimeDelta : nextTimeDecad);
}

long unsigned UpdateManager::computeNextEvent() const
{
	return _nextEvents + _howMany;
}

} /* namespace Kernel */
