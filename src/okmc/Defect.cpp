/*
 * Defect.cpp
 *
 *  Created on: Jan, 2012
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

#include "Defect.h"
#include "Interface.h"
#include "kernel/Mesh.h"
#include "kernel/MeshElement.h"
#include "kernel/MeshParam.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include <vector>

using std::vector;

namespace OKMC
{

//simplified to detect that materials are the same at both sides
//looking for the elements to be the same gives an error, because in diagonal movements the interface
//contains the adjacent element, not the diagonal.
//It assumes that the interface is either in the initial or the final box, so that it will not work for very tiny boxes...
Defect::INTERFACE_ACTIONS Defect::interactSurface(Kernel::SubDomain *pSub, Kernel::MeshElement *from, Kernel::MeshElement *to)
{
	vector<Particle *> dummy;

	const vector<Interface *> & interfacesF = from->getInterfaces();
	for(vector<Interface *>::const_iterator it=interfacesF.begin(); it!=interfacesF.end(); ++it)
	{
		if(((*it)->getElement(0)->getMaterial() == from->getMaterial() &&
				(*it)->getElement(1)->getMaterial() == to->getMaterial()) ||
		   ((*it)->getElement(1)->getMaterial() == from->getMaterial() &&
				   (*it)->getElement(0)->getMaterial() == to->getMaterial()))
		{
			if((*it)->canInteract(pSub, this))
			{
				(*it)->interact(pSub, this, 0, dummy);
				return INTERACTION_EXECUTED;
			}
			return INTERACTION_REJECTED;
		}
	}
	const vector<Interface *> & interfacesT = to->getInterfaces();
	for(vector<Interface *>::const_iterator it=interfacesT.begin(); it!=interfacesT.end(); ++it)
	{
		if(((*it)->getElement(0)->getMaterial() == from->getMaterial() &&
				(*it)->getElement(1)->getMaterial() == to->getMaterial()) ||
		   ((*it)->getElement(1)->getMaterial() == from->getMaterial() &&
				   (*it)->getElement(0)->getMaterial() == to->getMaterial()))
		{
			if((*it)->canInteract(pSub, this))
			{
				(*it)->interact(pSub, this, 0, dummy);
				return INTERACTION_EXECUTED;
			}
			return INTERACTION_REJECTED;
		}
	}
	return NO_INTERFACE;
}

//work around to obtain the particles in a const way for I/O
std::vector<const Particle *> Defect::getParticles() const
{
	Defect *This = const_cast<Defect *>(this);
	std::vector<const Particle *> parts;
	std::vector<Particle *> nc_parts = This->getParticles();
	for(std::vector<Particle *>::iterator it=nc_parts.begin(); it!=nc_parts.end(); ++it)
		parts.push_back(*it);
	return parts;
}

std::string Defect::name() const
{
	return Domains::global()->PM()->getIDName(getID());
}

//any particle present in list is to be removed from remove
void Defect::removeFromVector(std::vector<Particle *> &remove, const std::vector<Particle *> &list)
{
	for(std::vector<Particle *>::const_iterator it = list.begin(); it != list.end(); ++it)
		if(std::remove(remove.begin(), remove.end(), *it) != remove.end())
			remove.pop_back();
}

void Defect::breakPosition(Kernel::SubDomain *pSub, Kernel::Coordinates &c, Kernel::MeshElement *pEle, float breakDistance) const
{
	Kernel::Coordinates old(c);
	breakDistance -= pSub->_rng.rand()*_pDomain->_pMePar->_lambda[pEle->getMaterial()];
	float phi   = 2.*M_PI*pSub->_rng.rand();
	float theta =    M_PI*pSub->_rng.rand();
	c._x += sin(phi)*cos(theta)*breakDistance;
	c._y += sin(phi)*sin(theta)*breakDistance;
	c._z += cos(theta)*breakDistance;
	Kernel::MeshElement *pElement = pEle;
	if(pEle->getDomain()->_pMesh->checkMove(pElement ,c) != Kernel::Mesh::JUMP_OK)
		c = old;
}

void Defect::restart(std::ostream &os) const
{
	Event::restart(os);
}

}
