/*
 * CascadeEvent.h
 *
 *  Created on: May 6, 2015
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

#ifndef CASCADEEVENT_H_
#define CASCADEEVENT_H_

#include <fstream>
#include <vector>
#include "kernel/Event.h"
#include "kernel/Coordinates.h"

namespace Kernel { class Mesh; }

namespace Cascade {

class CascadeEvent: public Kernel::Event {
public:
	CascadeEvent(Kernel::Domain *pD, const std::vector<std::string> &f,
			bool bReact, bool periodic, bool voluminic, bool bDisplace, bool bCorrectX,
			const std::string &filename, float rate);
	virtual ~CascadeEvent();

	void      perform (Kernel::SubDomain *, unsigned);
	void      setIndex(unsigned, int idx) {_idx = idx; }
	unsigned  getNEvents() const { return 1; }
	int       getIndex(unsigned ) const { return _idx; }
	float     getRate(unsigned , float ) const { return _rate; }
	E_TYPE    getEType() const { return Event::CASCADE; }

private:
	int _idx;

	std::ifstream _pist;
	std::vector<int> _startPos; // vector storing all positions in the cascade file
	std::vector<int> _cPos;

	std::vector<std::string> _format;
	bool _bReact, _periodic, _voluminic, _bDisplace;
	bool _bCorrectX;
	float _rate;

	void create(Kernel::Mesh *,
			float x, float y, float z,
			const std::vector<std::string> &format, const std::vector<std::string> &fileLine,
			const Kernel::Coordinates &, const Kernel::Coordinates &,
			bool bReact, bool periodic, bool voluminic);
};

} /* namespace Cascade */

#endif /* CASCADEEVENT_H_ */
