/*
 * OKMCModel.cpp
 *
 *  Created on: May 11, 2011
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

#include "OKMCModel.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "kernel/Material.h"
#include "kernel/RateManager.h"
#include "kernel/Event.h"
#include "Interface.h"

namespace OKMC {

OKMCModel::OKMCModel(Kernel::Domain *pDomain)
{
	//Interfaces are inserted with the "left" element first. Left means negative sign
	LOWMSG("Interfacing...");
	Kernel::Mesh *pMesh = pDomain->_pMesh;
	for(Kernel::Mesh::iterator it=pMesh->begin(); it!=pMesh->end(); ++it)
	{
		Kernel::M_TYPE mt = it->getMaterial();
		unsigned ix, iy, iz;
		pMesh->getIndicesFromIndex(it->getIndex(), ix, iy, iz);
		if(ix)
		{
			Kernel::MeshElement *me2 = pMesh->getElement(pMesh->getIndexFromIndices(ix-1, iy, iz));
			if(me2->getMaterial() != mt)
			{
				Interface *pIF = new OKMC::Interface(pDomain, me2, &*it, 0);
				pDomain->_pRM->insert(pIF, me2);
			}
		}
		if(iy)
		{
			Kernel::MeshElement *me2 = pMesh->getElement(pMesh->getIndexFromIndices(ix, iy-1, iz));
			if(me2->getMaterial() != mt)
			{
				Interface *pIF = new OKMC::Interface(pDomain, me2, &*it, 1);
				pDomain->_pRM->insert(pIF, me2);
			}
		}
		if(iz)
		{
			Kernel::MeshElement *me2 = pMesh->getElement(pMesh->getIndexFromIndices(ix, iy, iz-1));
			if(me2->getMaterial() != mt)
			{
				Interface *pIF = new OKMC::Interface(pDomain, me2, &*it, 2);
				pDomain->_pRM->insert(pIF, me2);
			}
		}
	}
}

}

