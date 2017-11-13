/*
 * Particle.h
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

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "kernel/Coordinates.h"
#include "kernel/ParticleType.h"
#include <vector>

namespace Kernel
{
class SubDomain;
class Mesh;
class MeshElement;
}

namespace OKMC {

class Defect;

class Particle {
public:
	Particle(Kernel::P_TYPE type, const Kernel::Coordinates &c, Defect *, const Kernel::Coordinates &orig);
	Particle(std::istream &, Defect *);

	virtual ~Particle() {}

	Kernel::Coordinates getCoordinates() const { return _coord; }
	Kernel::Coordinates getOrig() const { return _orig; }
	void setOrig(const Kernel::Coordinates &o) { _orig = o; }
	void setCoordinates(const Kernel::Coordinates &c) { _coord = c; }
	Kernel::MeshElement * getElement() const { return _pElement; }
	Particle * getNext() const { return _next; }
	Defect * getDefect() const { return _pDefect; }
	void setDefect(Defect *p) { _pDefect = p; }
	Kernel::P_TYPE getPType() const { return _ptype; }
	void   setPType(Kernel::P_TYPE pt) { _ptype = pt; }

	void   selfdiffusion(Kernel::SubDomain *, Kernel::MeshElement *);
	void   restart(std::ostream &) const;

protected:
	Particle *_next, *_prev;
	Kernel::Coordinates _coord;
	Kernel::Coordinates _orig; //mainly used to compute diffusivities
	Kernel::MeshElement *_pElement;
	Defect *_pDefect;
	Kernel::P_TYPE _ptype;

    friend class Kernel::Mesh;
};

}

#endif /* PARTICLE_H_ */
