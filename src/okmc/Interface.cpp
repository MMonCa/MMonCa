/*
 * Interface.cpp
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

#include "Interface.h"
#include "kernel/MeshElement.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/Mesh.h"
#include "kernel/MeshParam.h"
#include "io/Diagnostic.h"
#include "MobileParticle.h"
#include "Cluster.h"
#include "kernel/RateManager.h"
#include "InterfaceParam.h"
#include "io/ParameterManager.h"
#include "domains/MCClient.h"
#include <cmath>

using std::vector;
using std::map;
using Kernel::M_TYPE;
using Kernel::P_TYPE;
using Kernel::POS_0;
using Kernel::MAX_PARTICLES;
using Kernel::V_TYPE;

namespace OKMC {

//m1 is the element with coordinates ix -1, iy -1 or iz -1, that is, with "negative" sign
Interface::Interface(Kernel::Domain *p, Kernel::MeshElement *m1, Kernel::MeshElement *m2, unsigned axis) : Defect(p), _axis(axis)
{
	_state=0;
	_me[0] = m1;
	_me[1] = m2;
	_myID._mt = (Domains::global()->PM()->isGasLike(m1->getMaterial()) ? m2->getMaterial() : m1->getMaterial());
	if(_pDomain->_pIFPar->_params[m1->getMaterial()][m2->getMaterial()] == 0 &&
	   _pDomain->_pIFPar->_params[m2->getMaterial()][m1->getMaterial()] == 0)
		ERRORMSG("Sorry, interface " << Domains::global()->PM()->getMaterialName(m1->getMaterial()) << " / " <<
				Domains::global()->PM()->getMaterialName(m2->getMaterial()) << " present, but not defined");
	int idx1 = m1->getIndex();
	int idx2 = m2->getIndex();
	unsigned ix1, iy1, iz1, ix2, iy2, iz2;
	_pDomain->_pMesh->getIndicesFromIndex(idx1, ix1, iy1, iz1);
	_pDomain->_pMesh->getIndicesFromIndex(idx2, ix2, iy2, iz2);
	_pDomain->_pMesh->getCorners(_me[1]->getIndex(), _coords_min, _coords_max);
	switch(axis)
	{
	case 0:
		assert(ix1 < ix2);
		_coords_max._x = _coords_min._x;
		_surfaceCm2 = (_coords_max._y - _coords_min._y)*(_coords_max._z - _coords_min._z)*1e-14;
		break;
	case 1:
		assert(iy1 < iy2);
		_coords_max._y = _coords_min._y;
		_surfaceCm2 = (_coords_max._x - _coords_min._x)*(_coords_max._z - _coords_min._z)*1e-14;
		break;
	case 2:
		assert(iz1 < iz2);
		_coords_max._z = _coords_min._z;
		_surfaceCm2 = (_coords_max._x - _coords_min._x)*(_coords_max._y - _coords_min._y)*1e-14;
		break;
	}
	int test = (ix1 != ix2? 1:0) + (iy1 != iy2? 1:0) + (iz1 != iz2? 1:0);
	if(test != 1)
		ERRORMSG("Interface internal error: Wrong elements!");
	m1->insert(this);
	m2->insert(this);
}

void Interface::restart(std::ostream &os) const
{
	os << _particles.size() << " " << _pDomain->_domain_number << " " << _me[0]->getIndex() << " " << _axis << " ";
	for(vector<Particle *>::const_iterator it = _particles.begin(); it != _particles.end(); ++it)
		(*it)->restart(os);
}

void Interface::restart(std::istream &is)
{
	unsigned size, domain, element, axis;
	is >> size >> domain >> element >> axis;
	// get correct interface
	Kernel::Domain *pDom = Domains::global()->getDomain(domain);
	Kernel::MeshElement *pME = pDom->_pMesh->getElement(element);
	Interface *pInt = 0;
	for(vector<Interface *>::iterator it = pME->getInterfaces().begin(); it!= pME->getInterfaces().end(); ++it)
		if((*it)->getAxis() == axis && (*it)->getElement(0) == pME)
			pInt = *it;
	if(!pInt)
		ERRORMSG("File corrupted when looking at interface definition!");
	for(unsigned i=0; i<size; ++i)
		pInt->putOnInterface(new Particle(is, pInt));
}

Interface::~Interface()
{
	_me[0]->remove(this);
	_me[1]->remove(this);

	for(vector<Particle *>::iterator it=_particles.begin(); it != _particles.end(); ++it)
	{
		_pDomain->_pMesh->remove(*it);
		delete *it;
	}
}

//emits at the "otherside"
void Interface::perform (Kernel::SubDomain *pSub, unsigned ev)
{
	if(ev < MAX_PARTICLES*2)
		emit(pSub, ev);
	else
		migrate(pSub, ev - MAX_PARTICLES*2);
}

void Interface::removeFromMap(P_TYPE family)
{
	map<P_TYPE, unsigned>::iterator it = _myID._pt.find(family);
	assert(it != _myID._pt.end());
	it->second--;
	_myID._pos.begin()->second--;
	Cluster::collapseMap(_myID);
}

//emit the particle in the side (both specified by ev)
void Interface::emit(Kernel::SubDomain *pSub, unsigned ev)
{
	unsigned side= ev % 2;
	unsigned otherside = side? 0:1;
	P_TYPE particle = ev/2;
	P_TYPE family   = Domains::global()->PM()->getFamily(particle);
	M_TYPE mt = _me[otherside]->getMaterial();
	float l = _pDomain->_pMePar->_lambda[mt];
	P_TYPE mt0 = Domains::global()->PM()->getMaterial(mt)._pt[0];
	P_TYPE mt1 = Domains::global()->PM()->getMaterial(mt)._pt[1];
	pSub->_evLog.performed(_me[otherside]->getMaterial(), Event::INTERFACE, mt, 0, particle, _state);

	Kernel::Coordinates coordi, origin;
	bool selfDefect = family == V_TYPE || family == mt0 || family == mt1;
	if(selfDefect) //Is, Vs or antisites
	{
		coordi = Kernel::Coordinates(
				_coords_min._x + (_coords_max._x - _coords_min._x)*pSub->_rng.rand(),
				_coords_min._y + (_coords_max._y - _coords_min._y)*pSub->_rng.rand(),
				_coords_min._z + (_coords_max._z - _coords_min._z)*pSub->_rng.rand());

		if(coordi._x >= _coords_max._x)
			coordi._x = _coords_max._x - 1e-3;
		if(coordi._y >= _coords_max._y)
			coordi._y = _coords_max._y - 1e-3;
		if(coordi._z >= _coords_max._z)
			coordi._z = _coords_max._z - 1e-3;
	}
	else
	{
		removeFromMap(family);
		for(vector<Particle *>::iterator itpart = _particles.begin(); itpart != _particles.end(); ++itpart)
			if((*itpart)->getPType() == family) //find particle and "emit" it
			{
				coordi = (*itpart)->getCoordinates();
				coordi[_axis] = _coords_min[_axis]; //starts from the pure interface
				origin = (*itpart)->getOrig();
				pSub->_pDomain->_pMesh->remove(*itpart);
				delete *itpart;
				*itpart = _particles.back();
				_particles.pop_back();
				break;
			}
		_pDomain->_pRM->update(this, _me[0]);
	}

	if(side)
		coordi[_axis] -= l * pSub->_rng.rand() + 0.001;
	else
		coordi[_axis] += l * pSub->_rng.rand() + 0.001;
	if(selfDefect)
		origin = coordi;
	MobileParticle *mp = new MobileParticle(pSub, particle, 0, _pDomain, coordi, origin); //self inserts
	if(mp->getElement() != _me[otherside])
		assert(false);
}

//perform one diffusion hop in the particle type specified by ev
void Interface::migrate(Kernel::SubDomain *pSub, unsigned ev)
{
	static int migMatrix[3][4] = { {1,1,2,2}, {0,0,2,2}, {0,0,1,1} };
	static int sigMatrix[4] = { 1, -1, 1, -1 };
	//pSub->_evLog.performed(_me[0]->getMaterial(), Event::INTERFACE, , 0, particle, _state);
	//pick a a random impurity of the correct type (ev);
	unsigned idx = unsigned(_myID._pt[ev] * pSub->_rng.rand());
	Particle *pPart = 0;
	unsigned count = 0;
	for(vector<Particle *>::iterator it=_particles.begin(); it!=_particles.end(); it++)
		if((*it)->getPType() == ev && count++ == idx)
		{
			pPart = *it;
			break;
		}
	assert(pPart);
	//chose lamdba
	M_TYPE mt0 = _me[0]->getMaterial();
	M_TYPE mt1 = _me[1]->getMaterial();
	const float l = _pDomain->_pIFPar->_params[mt0][mt1]->_lambda;
	const int axis = int(pSub->_rng.rand() * 4.);
	Kernel::Coordinates origDelta(0,0,0);
	origDelta[migMatrix[_axis][axis]] = l * sigMatrix[axis];
	Kernel::Coordinates c = pPart->getCoordinates();
	Kernel::Coordinates o = pPart->getOrig();
	Kernel::MeshElement *pEle = pPart->getElement();
	std::pair<P_TYPE, unsigned> stateID = std::make_pair(ev, 0);
	Kernel::Mesh::JUMP_ACTIONS action = _pDomain->_pMesh->jumpPosition(c, origDelta, pEle, &stateID, o, 1);
	if(action != Kernel::Mesh::JUMP_OK)
		return;
	//check that it goes into an interface
	vector<Interface *> & interfaces = pEle->getInterfaces();
	Interface *pnewInt = 0;
	for(vector<Interface *>::iterator it=interfaces.begin(); it!=interfaces.end(); ++it)
	{
		pnewInt = *it;
		if(pnewInt == this)
			break;
	}
	if(!pnewInt)
		return;
	pPart->setCoordinates(c);
	pPart->setOrig(o);
	//I have an interface!!
	if(pnewInt != this)
	{
		removeFromMap(ev);
		for(vector<Particle *>::iterator itpart = _particles.begin(); itpart != _particles.end(); ++itpart)
			if((*itpart)->getPType() == ev) //find particle and "emit" it
			{
				*itpart = _particles.back();
				_particles.pop_back();
				break;
			}
		_pDomain->_pMesh->remove(pPart);
		_pDomain->_pRM->update(this, _me[0]);
		pnewInt->insertParticle(pPart); //it inserts the particle in the mesh
		_pDomain->_pRM->update(pnewInt, pnewInt->_me[0]);
	}
}

float Interface::getRate(unsigned ev, float kT) const
{
	//emit
	if(ev < MAX_PARTICLES*2)
	{
		if(ev >= Domains::global()->PM()->getNParticles()*2)
			return 0;
		unsigned side=ev%2;
		unsigned otherside=side? 0:1;
		P_TYPE particle = ev/2;
		P_TYPE family = Domains::global()->PM()->getFamily(particle);
		M_TYPE mt0 = _me[side]->getMaterial();
		M_TYPE mt1 = _me[otherside]->getMaterial();
		P_TYPE pt0 = Domains::global()->PM()->getMaterial(mt1)._pt[0];
		P_TYPE pt1 = Domains::global()->PM()->getMaterial(mt1)._pt[1];
		bool selfDefect = family == V_TYPE || family == pt0 || family == pt1;               
		if(particle >= Domains::global()->PM()->getNParticles() || !Domains::global()->PM()->isParticleDefined(particle, mt1))
			return 0;
		if(_pDomain->_pIFPar->_params[mt0][mt1] == 0) //gasLike material
			return 0;
		float rate = _pDomain->_pIFPar->_params[mt0][mt1]->_trapMPProb[particle];                                
		rate *= _pDomain->_pIFPar->_params[mt0][mt1]->_arrEmitMP[particle](_me[otherside]).getRate(kT);
		if(selfDefect)
		{
			rate *= _surfaceCm2 * _pDomain->_pIFPar->_params[mt0][mt1]->_superSat[particle];
			return rate; //includes migration and barrier and recombination length
		}
		else
		{
			map<P_TYPE, unsigned>::const_iterator it= _myID._pt.find(family);
			if(it == _myID._pt.end())
				return 0;
			return rate * it->second; //includes migration and barrier and recombination length
		}
	}
	//migrate
	ev -= MAX_PARTICLES*2;
	map<P_TYPE, unsigned>::const_iterator it = _myID._pt.find(ev);
	M_TYPE mt0 = _me[0]->getMaterial();
	M_TYPE mt1 = _me[1]->getMaterial();
	if(_pDomain->_pIFPar->_params[mt0][mt1] == 0) //gasLike Material
		return 0;
	if(it != _myID._pt.end())
		return it->second * _pDomain->_pIFPar->_params[mt0][mt1]->_migration[ev].getRate(kT);
	return 0;
}

// insert a particle into the interface
// used when Mesh::
void Interface::insertParticle(Particle *part)
{
	const float epsilon = 0.01; //nm

	assert(part->getPType() < Domains::global()->PM()->getNFamilies());

	Kernel::Coordinates newC;
	Kernel::Coordinates oldC = static_cast<MobileParticle *>(part)->getCoordinates();
	newC._x = _coords_min._x + Domains::global()->client()->rand()*(_coords_max._x - _coords_min._x);
	newC._y = _coords_min._y + Domains::global()->client()->rand()*(_coords_max._y - _coords_min._y);
	newC._z = _coords_min._z + Domains::global()->client()->rand()*(_coords_max._z - _coords_min._z);

	if(_coords_max._x == _coords_min._x) newC._x = (oldC._x > _coords_max._x ? _coords_max._x + epsilon : _coords_max._x - epsilon);
	if(_coords_max._y == _coords_min._y) newC._y = (oldC._y > _coords_max._y ? _coords_max._y + epsilon : _coords_max._y - epsilon);
	if(_coords_max._z == _coords_min._z) newC._z = (oldC._z > _coords_max._z ? _coords_max._z + epsilon : _coords_max._z - epsilon);

	part->setCoordinates(newC);        // set new coordinates

	putOnInterface(part);
}

void Interface::putOnInterface(Particle *part)
{
	part->setDefect(this);             // set defect to current interface
	_particles.push_back(part);        // insert the particle in the interface
	_pDomain->_pMesh->insert(part, 0); // and in the mesh

	// finally, update interface properties
	P_TYPE family = part->getPType();
	if(_myID._pt.find(family) == _myID._pt.end())
		_myID._pt[family] = 0;
	_myID._pt[family]++;
	if(_myID._pos.empty())
		_myID._pos[POS_0] = 0;
	_myID._pos[POS_0]++;
}

// two possibilities:
// part is null. Then the defect interacts with the interface itself
// part is not null. The defect interacts with the part at the interface.
// Desorption probability is not used here (for Is and Vs) ,
// since it has already been taken into account in "Interface::canInteract".
// the parts vector needs to be accounted for!
Defect * Interface::interact(Kernel::SubDomain *pSub, Defect *dyn, Particle *sta, vector<Particle *> &parts)
{
	const float epsilon = 0.01; //nm
	unsigned def2 = 0;
	if(dyn->getEType() == MOBILEPARTICLE)
		def2 = static_cast<MobileParticle *>(dyn)->getPType();
	if(dyn->getEType() == CLUSTER)
		def2 = static_cast<Cluster *>(dyn)->getEDType();

	if (sta == 0) //interacts with the interface itself
	{
		if(dyn->getEType() == Event::MOBILEPARTICLE)
		{
			M_TYPE mt1 = dyn->getElement()->getMaterial();
			pSub->_reLog.reaction(mt1, INTERFACE, 0,_state, dyn->getEType(), def2, dyn->getState());
			M_TYPE mt0 = (mt1 == _me[1]->getMaterial()? _me[0]->getMaterial() : _me[1]->getMaterial());
			assert( mt1 == _me [1]->getMaterial() || mt1 == _me[0]->getMaterial());
			Kernel::Coordinates oldC = static_cast<MobileParticle *>(dyn)->getCoordinates();
			Kernel::Coordinates oldO = static_cast<MobileParticle *>(dyn)->getOrig();
			P_TYPE family = Domains::global()->PM()->getFamily(def2);
			removeFromVector(parts, dyn->getParticles());
			delete dyn;
			{
				P_TYPE pt0 = Domains::global()->PM()->getMaterial(mt1)._pt[0];
				P_TYPE pt1 = Domains::global()->PM()->getMaterial(mt1)._pt[1];
				bool selfDefect = family == V_TYPE || family == pt0 || family == pt1;
				if(selfDefect)
					return this;
			}
			//compute dopant concentration
			map<P_TYPE, unsigned>::iterator it=_myID._pt.find(family);
			float dopConc = (it == _myID._pt.end()? 0: it->second) / _surfaceCm2;
			float desorption = dopConc < _pDomain->_pIFPar->_params[mt0][mt1]->_desorptionThreshold[family] ?
					_pDomain->_pIFPar->_params[mt0][mt1]->_desorptionMPProbL[def2] :
					_pDomain->_pIFPar->_params[mt0][mt1]->_desorptionMPProbH[def2];
			if(pSub->_rng.rand() < desorption) //annihilation
				return this;
			assert(family < Domains::global()->PM()->getNFamilies());
			if(_coords_max._x == _coords_min._x) oldC._x = (oldC._x > _coords_max._x ? _coords_max._x + epsilon : _coords_max._x - epsilon);
			if(_coords_max._y == _coords_min._y) oldC._y = (oldC._y > _coords_max._y ? _coords_max._y + epsilon : _coords_max._y - epsilon);
			if(_coords_max._z == _coords_min._z) oldC._z = (oldC._z > _coords_max._z ? _coords_max._z + epsilon : _coords_max._z - epsilon);
			_particles.push_back(new Particle(family, oldC, this, oldO));
			if(it == _myID._pt.end())
				_myID._pt[family] = 0;
			_myID._pt[family]++;
			if(_myID._pos.empty())
				_myID._pos[POS_0] = 0;
			_myID._pos[POS_0]++;
			_pDomain->_pMesh->insert(_particles.back(), 0);
		}
		else if(dyn->getEType() == Event::CLUSTER)
		{
			M_TYPE mt1 = dyn->getElement()->getMaterial();
			M_TYPE mt0 = (mt1 == _me[1]->getMaterial()? _me[0]->getMaterial() : _me[1]->getMaterial());
			assert( mt1 == _me [1]->getMaterial() || mt1 == _me[0]->getMaterial());
			const Cluster *mc = static_cast<const Cluster *>(dyn);
			if(pSub->_rng.rand() >= _pDomain->_pIFPar->_params[mt0][mt1]->_desorptionMCProb[mc->getEDType()])
				return this;
			pSub->_reLog.reaction(mt1, INTERFACE, 0,_state, dyn->getEType(), def2, dyn->getState());
			pSub->_reLog.reactionInterface(mt1, def2, dyn->name());
			_pDomain->_pRM->remove(dyn, dyn->getElement());
			removeFromVector(parts, dyn->getParticles());
			delete dyn;
		}
		else
			ERRORMSG("Interaction of a non-implemented defect of type " << int(dyn->getEType()) << " with the interface!");
	}
	else
		WARNINGMSG("Interaction with surface particles not implemented yet");
	_pDomain->_pRM->update(this, _me[0]);
	return this;
}

//a particle interacts with a particle ON the interface. Returns false...
bool Interface::canInteract(Kernel::SubDomain *pSub, const Particle *dyn, const Particle *sta) const
{
	return false;
}

//a Defect interacts with the interface itself... it is OK
bool Interface::canInteract(Kernel::SubDomain *pSub, const Defect *def) const
{
	M_TYPE mt1 = def->getElement()->getMaterial();
	M_TYPE mt0 = (mt1 == _me[1]->getMaterial()? _me[0]->getMaterial() : _me[1]->getMaterial());
	assert( mt1 == _me[1]->getMaterial() || mt1 == _me[0]->getMaterial());

	if( def->getState() !=0)
		return false;
	switch(def->getEType())
	{
	case Event::MOBILEPARTICLE:
	{
		const MobileParticle *mp = static_cast<const MobileParticle *>(def);
		P_TYPE pt = mp->getPType();

		return _pDomain->_pIFPar->_params[mt0][mt1]->_interactWithMP[pt] &&    //if interaction is defined
				pSub->_rng.rand() < exp(-_pDomain->_pIFPar->_params[mt0][mt1]->_barrierMP[pt]/_pDomain->_pRM->getkT()) &&  //barrier term
				pSub->_rng.rand() < _pDomain->_pIFPar->_params[mt0][mt1]->_trapMPProb[pt];  // trap probability related to recombination length
	}
	case Event::CLUSTER:
	{
		const Cluster *mc = static_cast<const Cluster *>(def);
		bool it = _pDomain->_pIFPar->_params[mt0][mt1]->_interactWithMC[mc->getEDType()];
		float it2 = _pDomain->_pIFPar->_params[mt0][mt1]->_trapMCProb[mc->getEDType()];
		bool bTrap = pSub->_rng.rand() < it2;  //if recombination length is not defined, trapping is done by default
		return it && bTrap;
	}
	default:
		WARNINGMSG("IF: Interface interaction with defect " << Kernel::Event::getEName(def->getEType()) << " NOT supported");
		return false;
	}
}

void Interface::deletePart(Kernel::SubDomain *pSub, vector<Particle *> &parts)
{
	for(vector<Particle *>::iterator it=_particles.begin(); it != _particles.end(); ++it)
		if(*it == parts.back())
		{
			_pDomain->_pMesh->remove(*it);
			delete *it;
			*it = _particles.back();
			_particles.pop_back();
			removeFromMap(parts.back()->getPType());
			parts.pop_back();
			break;
		}
	_pDomain->_pRM->update(this, _me[0]);
}

}

