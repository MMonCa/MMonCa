/*
 * LatticeSite.cpp
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

#include "LatticeSite.h"
#include "kernel/Mesh.h"
#include "domains/Global.h"

namespace LKMC {

LatticeSite::LatticeSite(Kernel::M_TYPE basic, Kernel::P_TYPE type, const Kernel::Coordinates & c) : _coord(c), _type(type), _basicMat(basic)
{
	Domains::global()->getDomain(_coord)->_pMesh->insert(this, 0);
}

LatticeSite::LatticeSite(std::istream &is)
{
	is >> _coord >> _type >> _basicMat;
	Domains::global()->getDomain(_coord)->_pMesh->insert(this, 0);
}

void LatticeSite::restart(std::ostream &os) const
{
	os << _coord << " " << _type << " " << _basicMat << " ";
}

LatticeSite::~LatticeSite()
{
	Domains::global()->getDomain(_coord)->_pMesh->remove(this);
}

}
