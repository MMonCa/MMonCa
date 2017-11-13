/*
 * LatticeAtomBCC1.cpp
 *
 *  Created on: Dec 19, 2013
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

#include "LatticeAtomBCC1.h"
#include "LatticeBCC1Param.h"
#include "kernel/Domain.h"
#include "kernel/MeshElement.h"
#include "kernel/SubDomain.h"
#include "kernel/EventLog.h"
#include "kernel/RateManager.h"
#include "LKMCModel.h"

using Kernel::Coordinates;

namespace LKMC {

LatticeAtomBCC1::LatticeAtomBCC1(LASTATE per, Kernel::Domain *p, Kernel::M_TYPE basic, Kernel::P_TYPE type, const Coordinates &c) :
		LatticeAtom(per, p, basic, type, c)
{
	insertNeighbors();  //this function cannot be in the base class because it needs the virtual classes working
}

LatticeAtomBCC1::LatticeAtomBCC1(std::istream &is) : LatticeAtom(is)
{
	insertNeighbors();
}

void LatticeAtomBCC1::restart(std::ostream &os) const
{
	LatticeAtom::restart(os);
}

float LatticeAtomBCC1::getRate(unsigned eventType, float kT) const
{
	if(getPerformed() || eventType != 0)
		return 0;

	const LatticeBCC1Param * param = static_cast<const LatticeBCC1Param *>(_pDomain->_pLaPar[_basicMat]);

	int f = getCoordination(0);
	int s = getCoordination(1);
	if(f == 0 || f == 8)
		return 0;

	return param->_activation[f][s].getRate(kT);
}

void LatticeAtomBCC1::perform(Kernel::SubDomain *pSub, unsigned )
{
	pSub->_evLog.performed(0, Event::LATTICEATOM, LSDT_BCC, 0, _type, 0);
	_pDomain->_pRM->setDepthLA(_coord._x);
	_state = LS_PERFORMED;
	_pElement->updateLAStatus(pSub, MOD_Epitaxy);
	//update
	_pDomain->_pRM->update(this, _pElement);
	for(int i=0; i<THIRDN; ++i)
		if(_neighbors[i])
			_pDomain->_pRM->update(_neighbors[i], _neighbors[i]->getElement());
	updateME4Epitaxy(pSub);
}

void LatticeAtomBCC1::dumpXYZ(const std::vector<std::pair<int, Kernel::Coordinates> >&v, const std::string &name) const
{
	std::ofstream ll(name.c_str());

	ll << v.size() << std::endl << "dump\n";
	for(unsigned i=0; i<v.size(); ++i)
		ll << "W " << v[i].second._x*10 << " " << v[i].second._y*10 << " " << v[i].second._z*10 << std::endl;
	LOWMSG("File " << name << " created");
}

}
