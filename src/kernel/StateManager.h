/*
 * StateManager.h
 *
 *  Created on: Feb 20, 2013
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


#ifndef STATEMANAGER_H_
#define STATEMANAGER_H_

#include "Material.h"
#include "kernel/ParticleType.h"

namespace Kernel
{
class Domain;
class MeshElement;
class SubDomain;

class StateManager
{
public:
	StateManager(Kernel::Domain * p) : _pDomain(p) {};
	virtual ~StateManager() {}

	virtual unsigned changeState     (SubDomain *, P_TYPE, const MeshElement *) const = 0;
	virtual bool checkStateBarrier   (SubDomain *, P_TYPE, unsigned, const MeshElement *, const  MeshElement *) const = 0;
	virtual double bindingShift      (P_TYPE, P_POS, unsigned, const MeshElement *) const = 0;
	virtual bool     canInteract     (P_TYPE, unsigned, P_TYPE, unsigned, const MeshElement *) const = 0;
	virtual unsigned interactionState(P_TYPE, unsigned, P_TYPE, unsigned, const MeshElement *, P_TYPE) const = 0;
	virtual unsigned breakUpState    (P_TYPE, P_POS, unsigned, M_TYPE) const = 0;

protected:
	Kernel::Domain * _pDomain;
};
}


#endif /* STATEMANAGER_H_ */
