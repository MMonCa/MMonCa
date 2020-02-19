/*
 * EventLog.cpp
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

#include "EventLog.h"
#include "io/Diagnostic.h"
#include "Domain.h"
#include "domains/Splitter.h"
#include "io/ParameterManager.h"
#include "okmc/EDType.h"
#include "okmc/ClusterParam.h"
#include "lkmc/LatticeSite.h"
#include <iomanip>

using std::vector;
using std::map;
using std::stringstream;

namespace Kernel {

EventLog::EventLog()
{
	for(int i=0; i<Event::UNDEFINED_EVENT; ++i)
		_descriptions.push_back(std::vector<std::string>());
	_descriptions[Event::MOBILEPARTICLE].push_back("migrate");
	_descriptions[Event::MOBILEPARTICLE].push_back("break 0");
	_descriptions[Event::MOBILEPARTICLE].push_back("break 1");
	_descriptions[Event::MOBILEPARTICLE].push_back("emit I");
	_descriptions[Event::MOBILEPARTICLE].push_back("emit V");
	_descriptions[Event::MOBILEPARTICLE].push_back("state");
	_descriptions[Event::MOBILEPARTICLE].push_back("long hop");
	_descriptions[Event::MOBILEPARTICLE].push_back("rejected");

	_descriptions[Event::INTERFACE].push_back("Emission");
	_descriptions  [Event::CASCADE].push_back("Cascades");
	_descriptions  [Event::CELL].push_back("Create A");
	_descriptions  [Event::CELL].push_back("Erase A");	
	_descriptions  [Event::CELL].push_back("Create B");
	_descriptions  [Event::CELL].push_back("Erase B");	

	for(unsigned i=0; i<Domains::global()->getSplitter()->getLevels(); ++i)
	{
		std::stringstream ss;
		ss << i;
		_descriptions[Event::EMPTY].push_back(ss.str());
	}
}

void EventLog::performed(M_TYPE mt, Event::E_TYPE ev, unsigned defect_type, unsigned hash, unsigned defect_event, unsigned defect_state)
{
	while(_events[mt].size() <= unsigned(ev))
		_events[mt].push_back(vector<vector<vector<unsigned> > >());
	if(ev == Event::CLUSTER)
	{
		while(_clusters[mt].size() <= defect_type)
		{
			map<unsigned, vector<unsigned> > empty;
			_clusters[mt].push_back(empty);
		}
		map<unsigned, vector<unsigned> >::iterator it = _clusters[mt][defect_type].find(hash);
		if(it == _clusters[mt][defect_type].end())
		{
			_clusters[mt][defect_type][hash] = vector<unsigned>(MAX_PARTICLES+4,0);
			assert(_clusters[mt].size() > defect_type);
			assert(_clusters[mt][defect_type][hash].size() > defect_event);
			_clusters[mt][defect_type][hash][defect_event] = 1;
		}
		else
			(it->second)[defect_event]++;
	}
	else
	{
		while(_events[mt][ev].size() <= defect_type)
			_events[mt][ev].push_back(vector<vector<unsigned> >());
		while(_events[mt][ev][defect_type].size() <= defect_state)
			_events[mt][ev][defect_type].push_back(vector<unsigned>());
		while(_events[mt][ev][defect_type][defect_state].size() <= defect_event)
			_events[mt][ev][defect_type][defect_state].push_back(0);

		_events[mt][ev][defect_type][defect_state][defect_event]++;
	}
}

EventLog & EventLog::operator+=(const EventLog &log) //to be use when there are several domains.
{
	//_events[Iron][ExtendedDefect][I_TYPE][state][migration]
	//std::vector<std::vector<std::vector<std::vector<unsigned> > > >_events[MAX_MATERIALS];
	for(M_TYPE mt=0; mt<MAX_MATERIALS; ++mt)
	{
		while(_events[mt].size() < log._events[mt].size()) //defect not here.
			_events[mt].push_back(vector<vector<vector<unsigned> > >());

		for(unsigned ev=0; ev <  log._events[mt].size(); ++ev)
		{
			while(_events[mt][ev].size() < log._events[mt][ev].size()) // defect_type not here.
				_events[mt][ev].push_back(vector<vector<unsigned> >());
			for(unsigned defect_type=0; defect_type < log._events[mt][ev].size(); ++defect_type)
			{
				while(_events[mt][ev][defect_type].size() < log._events[mt][ev][defect_type].size()) // event type not here.
					_events[mt][ev][defect_type].push_back(vector<unsigned>());

				for(unsigned defect_state = 0; defect_state < log._events[mt][ev][defect_type].size(); ++defect_state)
				{
					while(_events[mt][ev][defect_type][defect_state].size() < log._events[mt][ev][defect_type][defect_state].size()) // defect state not here.
						_events[mt][ev][defect_type][defect_state].push_back(0);
					for(unsigned defect_event = 0; defect_event < log._events[mt][ev][defect_type][defect_state].size(); ++defect_event)
						_events[mt][ev][defect_type][defect_state][defect_event] += log._events[mt][ev][defect_type][defect_state][defect_event];
				}
			}
		}
	}
	//_clusters[material][311][hash][defect_state][event_type]
	//std::map<unsigned, std::vector<std::vector<unsigned> > > _clusters[MAX_MATERIALS+1];
	for(M_TYPE mt=0; mt<MAX_MATERIALS; ++mt)
	{
		while(_clusters[mt].size() <= log._clusters[mt].size())
		{
			map<unsigned, vector<unsigned> > empty;
			_clusters[mt].push_back(empty);
		}
		for(unsigned def_type = 0; def_type < log._clusters[mt].size(); ++def_type)
		{
			for(map<unsigned, vector<unsigned> >::const_iterator it=log._clusters[mt][def_type].begin();
					it != log._clusters[mt][def_type].end(); ++it) //for each hash in log
			{
				if(_clusters[mt][def_type].find(it->first) == _clusters[mt][def_type].end())
					_clusters[mt][def_type][it->first] = it->second;
				else
				{
					while(_clusters[mt][def_type][it->first].size() < it->second.size()) // defect_type not here.
						_clusters[mt][def_type][it->first].push_back(0);
					for(unsigned i = 0; i < _clusters[mt][def_type][it->first].size(); ++i)
						_clusters[mt][def_type][it->first][i] += it->second[i];
				}
			}
		}
	}
	return *this;
}

void EventLog::print() const
{
	LOWMSG("----------------- Event Log --------------");
	for(M_TYPE mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		for(unsigned ev = 0; ev < _events[mt].size(); ++ev)
		{
			bool bEvent = false;
			stringstream ssEv;
			ssEv << std::left << std::setw(18) << "Type";
			if(ev == Event::CLUSTER)
			{
				ssEv << std::right << std::setw(8) << "Mig";
				ssEv << std::right << std::setw(8) << "To";
				ssEv << std::right << std::setw(8) << "From";
				ssEv << std::right << std::setw(8) << "Rec";
				ssEv << std::left << std::setw(12) << " Emissions " << std::endl;
				for(unsigned def_type=0; def_type < _clusters[mt].size(); ++def_type)
					for(map<unsigned, vector<unsigned> >::const_iterator it=_clusters[mt][def_type].begin(); it!=_clusters[mt][def_type].end(); ++it)
					{
						stringstream clName;
						clName << Domains::global()->PM()->getEDName(mt, def_type) << "/";
						clName << Domains::global()->PM()->getIDName(Domains::global()->getDomain(0)->_pClPar->
								getParams(mt, def_type)->_hash.hash2cluster(it->first)) << " ";
						ssEv << std::left << std::setw(18) << clName.str();

						unsigned n=0;
						for(; n<3; ++n)
							if (it->second[MAX_PARTICLES+n])
							{
								bEvent = true;
								ssEv << std::right << std::setw(8)<< std::setprecision(2) << float(it->second[MAX_PARTICLES+n]);
							}
							else
								ssEv << std::right << std::setw(8) << " ";

						if(it->second[MAX_PARTICLES+n] != 0)
						{
							bEvent = true;
							ssEv << std::right << std::setw(8)<< std::setprecision(2) << float(it->second[MAX_PARTICLES+n]) << " ";
						}
						else
							ssEv << std::right << std::setw(9) << " ";

						for(unsigned pt=0; pt < MAX_PARTICLES; ++pt)
							if(it->second[pt] != 0)
							{
								bEvent = true;
								ssEv << it->second[pt] << '(' <<  Domains::global()->PM()->getParticleName(mt, pt) << ") ";
							}

						ssEv << std::endl;
					}
			}
			else if(ev == Event::LATTICEATOM)
			{
				stringstream ss;
				vector<std::string> descs[5];
				descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI1].push_back("Prec. 0");descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI1].push_back("Prec. 1");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI1].push_back("Prec. 2");descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI1].push_back("Migrat.");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI1].push_back("Etching");descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI1].push_back("Adsorpt");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI1].push_back("Desorpt");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI2].push_back("Prec. 0");descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI2].push_back("Prec. 1");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI2].push_back("Prec. 2");descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI2].push_back("Migrat.");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_EPI2].push_back("Etching");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_SPER].push_back("<100>_6");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_SPER].push_back("<100>_7");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_SPER].push_back("<100>_8");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_SPER].push_back("<100>_9");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_SPER].push_back("<100>_10");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_SPER].push_back("<110>");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_SPER].push_back("<111>");
				descs[LKMC::LatticeSite::LSDT_DIAMOND_SPER].push_back("none");
				for(unsigned defect_type=0; defect_type < _events[mt][ev].size(); ++defect_type)
					for(unsigned defect_state = 0; defect_state < _events[mt][ev][defect_type].size(); ++defect_state)
					{
						bool bPrint = false;
						stringstream ss;
						for(unsigned event_type = 0; event_type < descs[defect_type].size(); ++event_type)
							ss << std::right << std::setw(10) << descs[defect_type][event_type];
						ss << std::endl;
						ss << std::left << std::setw(18) << Event::getEName(Event::E_TYPE(ev));

						if(_events[mt][ev][defect_type].size())
						for(unsigned event_type = 0; event_type < _events[mt][ev][defect_type][defect_state].size(); ++event_type)
						{
							if(_events[mt][ev][defect_type][defect_state][event_type] != 0)
							{
								bPrint = true;
								if(ev == Event::INTERFACE)  //No state info given.
									ss << std::left << std::setw(4) << Domains::global()->PM()->getParticleName(mt, event_type) <<
									std::right << std::setw(10) << _events[mt][ev][defect_type][defect_state][event_type];
								else
									ss << std::right << std::setw(10) << _events[mt][ev][defect_type][defect_state][event_type];
							}
							else
							{
								if(ev != Event::INTERFACE)
									ss << std::right << std::setw(10) << " ";
								else
									ss << " ";
							}
						}
						if(bPrint)
						{
							bEvent = true;

							ssEv << ss.str() << std::endl;
						}
					}
			}
			else
			{
				for(unsigned event_type = 0; event_type < _descriptions[ev].size(); ++event_type)
					ssEv << std::right << std::setw(10) << _descriptions[ev][event_type];
				ssEv << std::endl;
				for(unsigned defect_type=0; defect_type < _events[mt][ev].size(); ++defect_type)
				{
					for(unsigned defect_state = 0; defect_state < _events[mt][ev][defect_type].size(); ++defect_state)
					{
						bool bId = true;
						bool bPrint = false;
						stringstream ss;
						if(ev == Event::MOBILEPARTICLE)
						{
							if(Domains::global()->PM()->getStates(mt,defect_type)>1)
								ss << std::left <<  Domains::global()->PM()->getParticleName(mt, defect_type);
							else
								ss << std::left << std::setw(18) << Domains::global()->PM()->getParticleName(mt, defect_type);
						}
						else if(ev == Event::INTERFACE)
						{
							ss << std::left << std::setw(18) << Domains::global()->PM()->getMaterialName(defect_type);
							bId = false;
						}
						else if(ev == Event::EMPTY)
						{
							ss << std::left << std::setw(18) << char('A' + defect_type);
							bId = false;
						}
						else
						{
							ss << std::left << std::setw(18) << Event::getEName(Event::E_TYPE(ev));
							bId = false;
						}

						if(Domains::global()->PM()->getStates(mt,defect_type) > 1 && bId &&
								Domains::global()->PM()->getStateName(mt, defect_type, defect_state) != "UNDEFINED_STATE")
						{
							std::string wholeName = Domains::global()->PM()->getParticleName(mt, defect_type) + "_" +  Domains::global()->PM()->getStateName(mt, defect_type, defect_state);
							ss<<  "_" <<Domains::global()->PM()->getStateName(mt, defect_type, defect_state)<<std::left << std::setw(18-wholeName.size())<< " ";
						}

						if(_events[mt][ev][defect_type].size())
							for(unsigned event_type = 0; event_type < _events[mt][ev][defect_type][defect_state].size(); ++event_type)
							{
								if(_events[mt][ev][defect_type][defect_state][event_type] != 0)
								{
									bPrint = true;
									if(ev == Event::INTERFACE)  //No state info given.
										ss << std::left << std::setw(4) << Domains::global()->PM()->getParticleName(mt, event_type) <<
										std::right << std::setw(10) << _events[mt][ev][defect_type][defect_state][event_type];
									else
										ss << std::right << std::setw(10) << _events[mt][ev][defect_type][defect_state][event_type];
								}
								else
								{
									if(ev != Event::INTERFACE)
										ss << std::right << std::setw(10) << " ";
									else
										ss << " ";
								}
							}
						if(bPrint)
						{
							bEvent = true;
							ssEv << ss.str() << std::endl;
						}
					}
				}
			}
			if(bEvent)
			{
				if(ev == Event::EMPTY)
					LOWMSG(" ----------------------------------------" <<
							std::left << std::setw(15) << " " <<
							std::right << std::setw(15) << Event::getEName(Event::E_TYPE(ev)));
				else
					LOWMSG(" ----------------------------------------" <<
							std::left << std::setw(15) << Domains::global()->PM()->getMaterialName(mt) <<
							std::right << std::setw(15) << Event::getEName(Event::E_TYPE(ev)));
				LOWMSG(ssEv.str());
			}
		}
}

}
