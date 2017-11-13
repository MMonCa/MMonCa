/*
 * ReactionLog.cpp
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

#include "ReactionLog.h"
#include "Domain.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "okmc/ClusterParam.h"
#include "okmc/EDType.h"
#include "kernel/Event.h"
#include <cassert>
#include <iomanip>

namespace Kernel {

using std::vector;
using std::map;
using std::stringstream;
using std::pair;
using std::string;

void ReactionLog::reaction(M_TYPE mt, Event::E_TYPE ev1, unsigned def1, unsigned st1, Event::E_TYPE ev2, unsigned def2, unsigned st2)
{
	// [mt][ev1]
	while(_reactions[mt].size() <= unsigned(ev1))
		_reactions[mt].push_back(vector<vector<vector<vector<vector<unsigned>  > > > >());
	// [mt][ev1][def1]
	while(_reactions[mt][ev1].size() <= def1)
		_reactions[mt][ev1].push_back(vector<vector<vector<vector<unsigned> > > >());
	// [mt][ev1][def1][st1]
	while(_reactions[mt][ev1][def1].size() <= st1)
		_reactions[mt][ev1][def1].push_back(vector<vector<vector<unsigned> > >());
	// [mt][ev1][def1][st1][ev2]
	while(_reactions[mt][ev1][def1][st1].size() <= unsigned(ev2))
		_reactions[mt][ev1][def1][st1].push_back(vector<vector<unsigned> >());
	// [mt][ev1][def1][st1][ev2][def2]
	while(_reactions[mt][ev1][def1][st1][ev2].size() <= def2)
		_reactions[mt][ev1][def1][st1][ev2].push_back(vector<unsigned>());
	// [mt][ev1][def1][st1][ev2][def2][st2]
	while(_reactions[mt][ev1][def1][st1][ev2][def2].size() <= st2)
		_reactions[mt][ev1][def1][st1][ev2][def2].push_back(0);
	_reactions[mt][ev1][def1][st1][ev2][def2][st2]++;
}

void  ReactionLog::reactionInterface(M_TYPE mt, unsigned defect1, const string &IDName)
{
	while(_reactionsInterface[mt].size() <= defect1)
		_reactionsInterface[mt].push_back(map<string, unsigned>());
	map<string, unsigned>::iterator it = _reactionsInterface[mt][defect1].find(IDName);
	if(it == _reactionsInterface[mt][defect1].end())
		_reactionsInterface[mt][defect1][IDName] = 1;
	else
		it->second++;
}

ReactionLog & ReactionLog::operator +=(const ReactionLog &log)
{
	//_reactions[mt][ev1][st1][def1][ev2][def2][st2]
	//        ev1        def1        st1         ev2         def2         st2
	// std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<unsigned> > > > > > _reactions[MAX_MATERIALS];

	for(M_TYPE mt=0; mt<MAX_MATERIALS; ++mt)
	{
		while(_reactions[mt].size() < log._reactions[mt].size()) //event type not here.
			_reactions[mt].push_back(vector<vector<vector<vector<vector<unsigned> > > > >());
		while(_reactionsInterface[mt].size() < log._reactionsInterface[mt].size())
			_reactionsInterface[mt].push_back(map<string, unsigned>());

		for(unsigned ev1=0; ev1 <  log._reactions[mt].size(); ++ev1)
		{
			while(_reactions[mt][ev1].size() < log._reactions[mt][ev1].size()) // defect_type not here.
				_reactions[mt][ev1].push_back(vector<vector<vector<vector<unsigned> > > >());
			for(unsigned def1=0; def1 < log._reactions[mt][ev1].size(); ++def1)
			{
				while(_reactions[mt][ev1][def1].size() < log._reactions[mt][ev1][def1].size()) // state not here.
					_reactions[mt][ev1][def1].push_back(vector<vector<vector<unsigned> > >());
				for (unsigned st1 = 0; st1 < log._reactions[mt][ev1][def1].size(); ++st1)
				{
					while(_reactions[mt][ev1][def1][st1].size() < log._reactions[mt][ev1][def1][st1].size()) // event type not here
						_reactions[mt][ev1][def1][st1].push_back(vector<vector<unsigned> >());
					for(unsigned ev2 = 0; ev2 < log._reactions[mt][ev1][def1][st1].size(); ++ev2)
					{
						while(_reactions[mt][ev1][def1][st1][ev2].size() < log._reactions[mt][ev1][def1][st1][ev2].size()) // defect_type not here.
							_reactions[mt][ev1][def1][st1][ev2].push_back(vector<unsigned>());
						for(unsigned def2=0; def2 < log._reactions[mt][ev1][def1][st1][ev2].size(); ++def2)
						{
							while(_reactions[mt][ev1][def1][st1][ev2][def2].size() < log._reactions[mt][ev1][def1][st1][ev2][def2].size()) // state not here.
								_reactions[mt][ev1][def1][st1][ev2][def2].push_back(0);
							for(unsigned st2 = 0; st2 <  log._reactions[mt][ev1][def1][st1][ev2][def2].size(); ++st2)
								_reactions[mt][ev1][def1][st1][ev2][def2][st2] += log._reactions[mt][ev1][def1][st1][ev2][def2][st2];
						}
					}
				}
			}
		}

		for(unsigned def1=0; def1 < log._reactionsInterface[mt].size(); ++def1)
		{
			for(map<string, unsigned>::const_iterator it = log._reactionsInterface[mt][def1].begin();
				it != log._reactionsInterface[mt][def1].end(); ++it)
				{
					map<string, unsigned>::iterator itMine = _reactionsInterface[mt][def1].find(it->first);
					if(itMine == _reactionsInterface[mt][def1].end())
						_reactionsInterface[mt][def1][it->first] = it->second;
					else
						itMine->second += it->second;
				}
		}
	}
	return *this;
}

void ReactionLog::print() const
{
	IO::ParameterManager *pPM = Domains::global()->PM();
	LOWMSG("----------------- Reaction Log --------------");
	for(M_TYPE mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		for(unsigned ev1 = 0; ev1 < _reactions[mt].size(); ++ev1)
		{
			bool bEvent = false;
			stringstream ssEv;
			for(unsigned def1=0; def1 < _reactions[mt][ev1].size(); ++def1)
			{
				bool bPrint = false;
				stringstream ss;
				ss << std::left << std::setw(10) << Event::getEName(Event::E_TYPE(ev1));
				//so far it has written something like
				// MobileParticle
				for(unsigned st1 = 0; st1 < _reactions[mt][ev1][def1].size(); ++st1)
					for(unsigned ev2 = 0; ev2 < _reactions[mt][ev1][def1][st1].size(); ++ev2)
						for(unsigned def2 = 0; def2 < _reactions[mt][ev1][def1][st1][ev2].size(); ++def2)
							for(unsigned st2 = 0; st2 < _reactions[mt][ev1][def1][st1][ev2][def2].size(); ++st2)
								if(_reactions[mt][ev1][def1][st1][ev2][def2][st2] != 0)  //finally there is something
								{
									bPrint = true;
									string state1 = "", state2 = "";
									if(ev1 == Event::MOBILEPARTICLE && pPM->getStates(mt, def1)!=1)
										state1 = "_"+ pPM->getStateName(mt, def1, st1);
									if(ev1 == Event::MOBILEPARTICLE && pPM->getStates(mt, def2)!=1)
										state2 = "_" + pPM->getStateName(mt, def2, st2);

									string defName = pPM->getDefName(mt, ev1, def1) + state1 + "+" +
											Domains::global()->PM()->getDefName(mt, ev2, def2) + state2;
									// has ICluster+V in def1Name
									ss << " " << std::left << std::setw(15) << defName <<
											std::right << std::setw(8) << _reactions[mt][ev1][def1][st1][ev2][def2][st2];
								}
				if(bPrint)
				{
					bEvent = true;
					ssEv << ss.str() << std::endl;
				}
			}
			if(bEvent)
			{
				LOWMSG(" ----------------------------------------" <<
						std::left << std::setw(15) << Domains::global()->PM()->getMaterialName(mt) <<
						std::right << std::setw(15) << Event::getEName(Event::E_TYPE(ev1)));
				LOWMSG(ssEv.str());
			}
		}
}

void ReactionLog::printInterface() const
{
	IO::ParameterManager *pPM = Domains::global()->PM();
	LOWMSG("----------------- Interface reaction Log --------------");
	for(M_TYPE mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		for(unsigned def1 = 0; def1 < _reactionsInterface[mt].size(); ++def1)
		{
			stringstream ss;
			for(std::map<string, unsigned>::const_iterator it=_reactionsInterface[mt][def1].begin();
					it != _reactionsInterface[mt][def1].end(); ++it)
				ss << std::left  << std::setw(10) << pPM->getDefName(mt, Event::CLUSTER, def1) << " "
				   << std::left  << std::setw(10) << it->first << " "
				   << std::right << std::setw(10) << it->second << std::endl;
			if(ss.str().size())
				LOWMSG(ss.str());
		}
}

}

