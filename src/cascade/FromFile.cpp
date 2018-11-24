/*
 * FromFile.cpp
 *
 *  Created on: Jul 21, 2011
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

#include "FromFile.h"
#include "CascadeEvent.h"
#include "io/Diagnostic.h"
#include "domains/MCClient.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "kernel/RateManager.h"
#include <fstream>

using std::string;
using std::vector;

namespace Cascade {

FromFile::FromFile()
{
	Domains::global()->client()->beginInsert();
}

FromFile::~FromFile()
{
	Domains::global()->client()->endInsert();
}

void FromFile::operator()(const string &fileName,
		vector<string> &format,
		float fluence,	float flux,
		string defects,   //procedure to obtain the defect name for clusters. (V4 is ICluster, 111?)
		bool bDisplace,   //do not move y and z coordinate
		bool periodic, bool bReact,
		float tempC,      //-1 means previous temperature
		bool voluminic,   // false means that cascades are just in the surface
		bool bCorrectX)   //set the init of an interface as x=0;
{
	std::ifstream pist(fileName.c_str());
	if(pist.fail())
		ERRORMSG("cascade: file " << fileName << " not found");

	if(fluence == 0)
		return;

	std::vector<CascadeEvent *> events;
	for(unsigned domain=0; domain < Domains::global()->getDomains(); ++domain)
	{
		Kernel::Domain * pDomain = Domains::global()->getDomain(domain);
		Kernel::Coordinates m,M;
		pDomain->_pMesh->getDomain(m,M);
		Kernel::MeshElement *pME = pDomain->_pMesh->getElement(0);
		float const howMany = fluence*1e-14 *(M._y - m._y)*(M._z - m._z);
		float rate = flux*howMany / fluence;
		events.push_back(new CascadeEvent(pDomain, format, bReact, periodic, voluminic, bDisplace, bCorrectX, fileName, rate));
		pDomain->_pRM->insert(events.back(), pME); //in the first element.
		if(flux == 0) //instantaneus...
		  {
		    unsigned hM = unsigned(howMany + pDomain->_rng_dom.rand());
		    for(unsigned i=0; i<hM; ++i)
		      events.back()->perform(pDomain->_pRM->getSubDomain(pME->getSubDomainIdx()), 0);
		  }
	}

	if(tempC != -1)
		Domains::global()->setTempK(tempC+273.15);
	if(flux != 0)
		Domains::global()->anneal(fluence/flux, false, false, 0);
	for(unsigned domain=0; domain < Domains::global()->getDomains(); ++domain)
	{
		Kernel::Domain * pDomain = Domains::global()->getDomain(domain);
		Kernel::MeshElement *pME = pDomain->_pMesh->getElement(0);
		pDomain->_pRM->remove(events[domain], pME);
		delete events[domain];
	}
}

}
