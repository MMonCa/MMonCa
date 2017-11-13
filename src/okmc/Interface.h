/*
 * Interface.h
 *
 *  Created on: Feb 24, 2011
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

#ifndef INTERFACE_H_
#define INTERFACE_H_

#include "Defect.h"
#include "kernel/Coordinates.h"
#include <vector>

namespace Kernel { class MeshElement; class Domain; }

namespace OKMC {

class Particle;

class Interface: public Defect
{
public:
	Interface(Kernel::Domain *, Kernel::MeshElement *, Kernel::MeshElement *, unsigned axis);
	virtual ~Interface();

	unsigned getAxis() const { return _axis; }
	Kernel::MeshElement * getElement(unsigned i) const { return _me[i]; }
	//it is an event
	virtual void perform (Kernel::SubDomain *, unsigned eventType);
    virtual void setIndex(unsigned eventType, int idx) { _idx[eventType] = idx; }
    virtual unsigned getNEvents() const { return Kernel::MAX_PARTICLES*2+Kernel::MAX_IMPURITIES; }
    virtual int  getIndex(unsigned eventType) const    { return _idx[eventType];}
    virtual float getRate(unsigned eventType, float kT) const;
    virtual E_TYPE getEType() const { return INTERFACE; }

	void insertParticle(Particle *);

    virtual Defect * interact(Kernel::SubDomain *, Defect *dyn, Particle *sta, std::vector<Particle *> &);
    virtual bool canInteract (Kernel::SubDomain *, const Particle *dyn, const Particle *sta) const;
   	virtual Kernel::ID getID() const  { return _myID; }
   	virtual std::vector<Particle *> getParticles() { return _particles; }
   	virtual Kernel::MeshElement * getElement() const { return _me[0]; }

   	virtual unsigned getState() const {return 0;}
   	virtual void deletePart(Kernel::SubDomain *, std::vector<Particle *> &);

    bool canInteract (Kernel::SubDomain *, const Defect *dyn) const; //defect with interface, only for interfaces.

    virtual void     restart(std::ostream &) const;
    static  void     restart(std::istream &);

private:
   	void emit   (Kernel::SubDomain *, unsigned ev);
   	void migrate(Kernel::SubDomain *, unsigned ev);
   	void removeFromMap(Kernel::P_TYPE family);
   	void putOnInterface(Particle *);

	unsigned _axis;
	Kernel::MeshElement *_me[2];
	unsigned _idx[Kernel::MAX_PARTICLES*2+Kernel::MAX_IMPURITIES]; //emission and migration
	Kernel::Coordinates _coords_min, _coords_max;
	float _surfaceCm2;
	std::vector<Particle *> _particles;
	Kernel::ID _myID;

	friend class Kernel::MeshElement;
};

}

#endif /* INTERFACE_H_ */
