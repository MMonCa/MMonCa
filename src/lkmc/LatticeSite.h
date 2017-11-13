/*
 * LatticeSite.h
 *
 *  Created on: May 30, 2011
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

#ifndef LATTICESITE_H_
#define LATTICESITE_H_

#include "kernel/Coordinates.h"
#include "kernel/ParticleType.h"

namespace Kernel { class Mesh; class MeshElement; class Domain; }

namespace LKMC {

class LatticeSite {
public:
	LatticeSite(Kernel::M_TYPE , Kernel::P_TYPE , const Kernel::Coordinates &);
	LatticeSite(std::istream &);
	virtual ~LatticeSite();

	Kernel::MeshElement * getElement() const { return _pElement; }
	LatticeSite * getNext() const { return _next; }
	Kernel::Coordinates getCoordinates() const { return _coord; }
	Kernel::P_TYPE getPType() const { return _type; }
	Kernel::M_TYPE getBasicMat() const { return _basicMat; }
	void setPType(Kernel::P_TYPE pt) { _type = pt; }
	void restart(std::ostream &) const;

	enum { LSDT_BCC, LSDT_FCC, LSDT_DIAMOND_SPER, LSDT_DIAMOND_EPI1, LSDT_DIAMOND_EPI2 };

protected:
	Kernel::Coordinates _coord;
	LatticeSite *_next, *_prev;
	Kernel::MeshElement * _pElement;
	Kernel::P_TYPE _type;
	Kernel::M_TYPE _basicMat;

	friend class Kernel::Mesh;
};

}

#endif /* LATTICESITE_H_ */
