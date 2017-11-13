/*
 * ReactionLog.h
 *
 *  Created on: Jul 21, 2011
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

#ifndef REACTIONLOG_H_
#define REACTIONLOG_H_

#include "Event.h"
#include "kernel/ParticleType.h"
#include "Material.h"
#include <vector>
#include <string>
#include <map>

namespace Kernel {

class ReactionLog {
public:
	ReactionLog() {}
	~ReactionLog() {}

	void print() const;
	void printInterface() const;
	void reaction(M_TYPE, Event::E_TYPE, unsigned defect1, unsigned state1, Event::E_TYPE, unsigned defect2, unsigned state2);
	void reactionInterface(M_TYPE, unsigned defect1, const std::string &IDName);
	ReactionLog & operator+=(const ReactionLog &);

private:
	//_reactions[mt][ev1][def1][st1][ev2][def2][st2]
	//     ev1        def1        st1         ev2         def2        st2
	std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<unsigned> > > > > > _reactions[MAX_MATERIALS];
	//_reactionsInteface[mt][def][IDName]  For clusters + interface only!
	std::vector<std::map<std::string, unsigned> > _reactionsInterface[MAX_MATERIALS];
};

}

#endif /* REACTIONLOG_H_ */
