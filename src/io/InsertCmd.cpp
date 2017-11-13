/*
 * InsertCmd.cpp
 *
 *  Created on: May 19, 2011
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

#include "InsertCmd.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/RateManager.h"
#include "kernel/Mesh.h"
#include "domains/MCClient.h"
#include "io/ParameterManager.h"
#include "okmc/MobileParticle.h"
#include "okmc/EDTypeIrregular.h"
#include "okmc/Cluster.h"
#include "domains/Global.h"

using std::string;
using Kernel::Coordinates;
using std::vector;

namespace IO {

// InsertCmd::InsertCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv, false)
InsertCmd::InsertCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv)
{
}

int InsertCmd::operator()()
{
	Coordinates c(getCoordinates("coord"));
	Kernel::Domain * pDomain = Domains::global()->getDomain(c);
	if(specified("interface"))
	{
		string ptName = getString("particle");
		Domains::global()->client()->createIF(ptName, c);
		return TCL_OK;
	}
	if(specified("particle"))
	{
		bool bReact = !specified("do.not.react");
		string ptName = getString("particle");
		string state  = specified("state")	?	getString("state") : "";
		Domains::global()->client()->createMP(ptName,state,c, bReact);
		return TCL_OK;
	}
	//extended defect
	if(specified("defect"))
	{
		string defName = getString("defect");
		bool bReact = !specified("do.not.react");
		Kernel::M_TYPE mt = pDomain->_pMesh->getElement(pDomain->_pMesh->getIndexFromCoordinates(c))->getMaterial();
		unsigned defNumber = pDomain->_pClPar->getDefectNumber(mt, defName);
		if(defNumber != pDomain->_pClPar->defectSize(mt))  //extended defect
		{
			Kernel::MeshElement *pEle = pDomain->_pMesh->getElement(pDomain->_pMesh->getIndexFromCoordinates(c));
			string ID = getString("ID");
			Kernel::ID theMap = Domains::global()->PM()->getID(pEle->getMaterial(), ID);
			if(theMap._pt.size() == 0)
				ERRORMSG("insert defect: Cannot understand ID " << ID);
			//check that we have the params needed for the defect
			if(pDomain->_pClPar->getParams(defNumber, theMap))
			{
				Kernel::SubDomain *pSub = pDomain->_pRM->getSubDomain(pEle->getSubDomainIdx());
				OKMC::Cluster *pCl = new OKMC::Cluster(pSub, pDomain, defNumber, c);
				pCl = pCl->growDefect(pSub, theMap);
				if(pCl)
				{
					pDomain->_pRM->insert(pCl, pCl->getElement());
					if(bReact)
					{
						std::vector<OKMC::Particle *> reacting, constituents = pCl->getParticles();
						for(std::vector<OKMC::Particle *>::iterator it=constituents.begin(); it!=constituents.end(); ++it)
							pDomain->_pMesh->getInteractions(pSub, *it, reacting);
						if(reacting.size())
						{
							OKMC::Particle *thePart = reacting[pSub->_rng.rand()*reacting.size()];
							std::vector<OKMC::Particle *> dummy;
							thePart->getDefect()->interact(pSub, pCl, thePart, dummy);
						}
					}
				}
				return TCL_OK;
			}
			else
				ERRORMSG("Cluster to be inserted is not defined in the parameters!");
			return TCL_OK;
		}
		else
			ERRORMSG("Defect name " << defName << " is not supported");
	}
	Tcl_AppendResult(_pTcl, "Option is missing: particle, defect", 0);
	return TCL_ERROR;
}

}
