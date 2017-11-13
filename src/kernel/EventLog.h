/*
 * EventLog.h
 *
 *  Created on: Apr 11, 2011
 *
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

#ifndef EVENTLOG_H_
#define EVENTLOG_H_

#include "Event.h"
#include "kernel/ParticleType.h"
#include "Material.h"
#include <vector>
#include <string>
#include <map>

namespace Kernel {

class EventLog {
public:
	EventLog();
	~EventLog() {};

	void print() const;
	void performed(M_TYPE mt, Event::E_TYPE ev, unsigned defect_type, unsigned hash, unsigned defect_event, unsigned defect_state);

	EventLog & operator+=(const EventLog &); //to be use when there are several domains.

private:
    //mt multicluster defect_type (for mobile parts) state event_type(mig)
	std::vector<std::vector<std::vector<std::vector<unsigned> > > > _events[MAX_MATERIALS];
	//defect_type        hash      defect_event
	std::vector<std::map<unsigned, std::vector<unsigned> > >_clusters[MAX_MATERIALS+1];
	std::vector<std::vector<std::string> > _descriptions;
};

}

#endif /* EVENTLOG_H_ */
