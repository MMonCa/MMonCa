/*
 * Defect.h
 *
 *  Created on: Mar 1, 2011
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

#ifndef DEFECT_H_
#define DEFECT_H_

#include "kernel/Event.h"
#include "kernel/Coordinates.h"
#include "kernel/ParticleType.h"
#include <string>
#include <vector>
#include <map>

namespace Kernel { class Domain; class MeshElement; }

namespace OKMC
{
class Particle;

class Defect: public Kernel::Event
{
public:
	Defect(Kernel::Domain *p) : Event(p) {}
	Defect(std::istream &is) : Event(is) {}

	virtual ~Defect() {}

	//sta->getDefect is the same as the defect being called.
	//the dynamic defect (moves and) interacts with the defect in the static particle. Also, some particles might be destroyed.
	// All the particles accounted in the interaction are to be removed from the vector.
	virtual Defect * interact(Kernel::SubDomain *, Defect *dyn, Particle *sta, std::vector<Particle *> &) = 0;
	// Check whether the interaction is possible including the capture radius!
	virtual bool canInteract(Kernel::SubDomain *, const Particle *dyn, const Particle *sta) const = 0;
	virtual Kernel::ID getID() const = 0; //slow, use just for interfacing.
	virtual std::vector<Particle *> getParticles() = 0; //for I/O operations
	std::vector<const Particle *> getParticles() const;
	virtual Kernel::MeshElement * getElement() const = 0;
	virtual void deletePart(Kernel::SubDomain *pSub, std::vector<Particle *> &) = 0; //deletes the back particle of the vector and takes care
		// of updating any of the other included particles with the new pointers if touched or deleted.
	enum INTERFACE_ACTIONS { NO_INTERFACE, INTERACTION_REJECTED, INTERACTION_EXECUTED };

	virtual unsigned getState() const=0;
	std::string name() const;

    virtual void     restart(std::ostream &) const;

protected:
	void removeFromVector(std::vector<Particle *> &remove, const std::vector<Particle *> &list);
	INTERFACE_ACTIONS interactSurface(Kernel::SubDomain *, Kernel::MeshElement *from, Kernel::MeshElement *to);  //this will react with a proper interface
	void  breakPosition(Kernel::SubDomain *, Kernel::Coordinates &, Kernel::MeshElement *, float breakDistance) const;
	unsigned _state;
};

}

#endif /* DEFECT_H_ */
