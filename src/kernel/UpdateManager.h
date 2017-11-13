/*
 * UpdateManager.h
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

#ifndef UPDATEMANAGER_H_
#define UPDATEMANAGER_H_

#include <string>

namespace Kernel {

class UpdateManager {
public:
	UpdateManager(const std::string &);
	bool operator()(double time, long unsigned events, double depth) const;
	double getNextTime() const { return _nextTime; }
	long unsigned getNextEvents() const { return _nextEvents; }
	void print(bool b) { _bPrint = b; }

private:
	double computeNextTime() const;
	long unsigned computeNextEvent() const;
	long unsigned _howMany;
	double _deltaT;
	double _decadesT;
	double _minT;
	const std::string _tag;
	bool _bPrint;

	mutable long unsigned _oldEvents;
	mutable long unsigned _nextEvents;
	mutable double _nextTime;
};

} /* namespace Kernel */
#endif /* UPDATEMANAGER_H_ */
