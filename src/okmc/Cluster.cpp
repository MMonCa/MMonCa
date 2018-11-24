/*
 * Cluster.cpp
 *
 *  Created on: May 30, 2011
 *  Modified on: August 2013
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

#include "Cluster.h"
#include "ClusterReactionParam.h"
#include "Particle.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/MeshElement.h"
#include "kernel/Mesh.h"
#include "MobileParticle.h"
#include "MobileParticleParam.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "kernel/RateManager.h"
#include "kernel/MeshParam.h"
#include "EDType.h"
#include "Interface.h"
#include <algorithm>

using Kernel::Domain;
using std::map;
using std::vector;
using std::string;
using Kernel::Coordinates;
using Kernel::P_POS;
using Kernel::POS_0;
using Kernel::POS_1;
using Kernel::POS_I;
using Kernel::POS_V;
using Kernel::P_TYPE;
using Kernel::V_TYPE;
using Kernel::ID;
using Kernel::UNDEFINED_TYPE;
using Kernel::MAX_PARTICLES;

namespace OKMC {

Cluster::Cluster(Kernel::SubDomain *pSub, Kernel::Domain *pDom, unsigned edType, MobileParticle *p1, MobileParticle *p2) : Defect(pDom)
{
	_edtype = edType;
	_myParams = 0;
	_state = 0;
	_map._mt = p1->getElement()->getMaterial();
	_center = p1->getCoordinates();
	_pElement = p1->getElement();
	setAxes(pSub); //needs mt and _pElement to be set
	_radius = 0;
	addParticle(pSub, p1);
	fuzzModel(pSub, p1->getPType(), p1->getCoordinates(), false); //cannot erase the cluster
	addParticle(pSub, p2);
	fuzzModel(pSub, p2->getPType(), p2->getCoordinates(), false);
	_myParams = getMyParam();
	delete p1;
	delete p2;
}

Cluster::Cluster(Kernel::SubDomain *pSub, Kernel::Domain *pDom, unsigned edType, const Coordinates &coord) : Defect(pDom)
{
	_edtype = edType;
	_state = 0;
	_center = coord;
	_radius = 0;
	_pElement = _pDomain->_pMesh->getElement(_pDomain->_pMesh->getIndexFromCoordinates(_center));
	_map._mt = _pElement->getMaterial();
	setAxes(pSub); //need mt to be set
}

Cluster::Cluster(std::istream &is) : Defect(is)
{
	is >> _edtype >> _state >> _center;
	_pElement = _pDomain->_pMesh->getElement(_pDomain->_pMesh->getIndexFromCoordinates(_center));
	_map._mt = _pElement->getMaterial();
	_radius = 0;
	setAxes(_pDomain->_pRM->getSubDomain(0)); //anyone should work

	//map
	unsigned nPos, nPt;
	is >> nPos >> nPt;
	for(unsigned i=0; i<nPos; ++i)
	{
		unsigned pos, n;
		is >> pos >> n;
		_map._pos[P_POS(pos)] = n;
	}
	for(unsigned i=0; i<nPt; ++i)
	{
		unsigned pt, n;
		is >> pt >> n;
		_map._pt[pt] = n;
	}
	unsigned nParticles;
	is >> nParticles;
	for(unsigned i=0; i<nParticles; ++i)
	{
		Particle *pPart = new Particle(is, this);
		_particles.push_back(pPart);
		updateRadius(pPart, +1);
		_pDomain->_pMesh->insert(pPart, 0);
	}
	_myParams = getMyParam();
}

void Cluster::restart(std::ostream &os) const
{
	Defect::restart(os);

	os << _edtype << " " << _state << " " << _center << " ";
	//map
	os << _map._pos.size() << " " << _map._pt.size() << " ";
	for(map<P_POS, unsigned>::const_iterator it=_map._pos.begin(); it != _map._pos.end(); ++it)
		os << int(it->first) << " " << it->second << " ";
	for(map<P_TYPE, unsigned>::const_iterator it=_map._pt.begin(); it != _map._pt.end(); ++it)
		os << int(it->first) << " " << it->second << " ";
	//particles
	os << _particles.size() << " ";
	for(unsigned i=0; i<_particles.size(); ++i)
		_particles[i]->restart(os);
}

/*
 * Creates a vector of particle names, picks up one randomly and uses surface to create it.
 * If it cannot be created, creates another cluster on top of it with the reminder.
 * if bAbort is true, an incomplete defect will not be created.
 *
 * Should NOT be included in RM when calling it!
 *
 */
Cluster * Cluster::growDefect(Kernel::SubDomain *pSub, const ID &theMap)
{
	ID tempMap = theMap;
	vector<P_TYPE> toInsert;
	vector<P_POS> toInsertPos;
	for(map<P_TYPE, unsigned>::const_iterator it=theMap._pt.begin(); it!=theMap._pt.end(); ++it)
		for(unsigned i=0; i<it->second; ++i)
			toInsert.push_back(it->first);
	for(map<P_POS, unsigned>::const_iterator it=theMap._pos.begin(); it!=theMap._pos.end(); ++it)
			for(unsigned i=0; i<it->second; ++i)
				toInsertPos.push_back(it->first);
	if(toInsertPos.size() != toInsert.size())
		ERRORMSG("Wrong ID in Cluster::growDefect.");
	vector<std::pair<Kernel::Coordinates, P_TYPE> > fuzzArray;
	while(toInsert.size())
	{
		unsigned idx = unsigned(pSub->_rng.rand()*toInsert.size());
		P_TYPE pt = toInsert[idx];
		P_POS pos = toInsertPos[idx];
		P_TYPE part = Domains::global()->PM()->getParticleForCluster(_map._mt, pt, pos);
		Particle *pP = new Particle(part, _center, this, _center);
		fuzzArray.push_back(std::make_pair(addParticle(pSub, pP), part));
		delete pP;
		tempMap._pt[pt]--;
		tempMap._pos[pos]--;
		toInsert.erase(toInsert.begin()+idx);
		toInsertPos.erase(toInsertPos.begin()+idx);
	}
	vector<Particle *> dummy;
	Defect * pDef = checkSize(pSub, false, dummy);
	if(pDef)
	{
		_pDomain->_pRM->insert(this, _pElement);
		for(vector<std::pair<Kernel::Coordinates, P_TYPE> >::iterator iter = fuzzArray.begin(); iter != fuzzArray.end(); ++iter)
			if(fuzzModel(pSub, iter->second, iter->first, true))
				return 0;
		_pDomain->_pRM->remove(this, _pElement);
	}
	return (pDef == this? this : 0);
}

void Cluster::setAxes(Kernel::SubDomain *pSub)
{
	const EDType *edt = _pDomain->_pClPar->getParams(_map._mt, _edtype);
	do
	{
		Coordinates axes[3];
		//changes the sign
		for(int j=0; j<3; ++j) //x, y or z
		{
			int sign = (pSub->_rng.rand() < 0.5? 1:-1);
			for(int i=0; i<3; ++i)
				axes[i][j] = sign*edt->_axes[i][j];
		}
		//permutation. There are six possible ones
		int permutations = int(pSub->_rng.rand()*6.)+1;
		vector<unsigned> perm;
		perm.push_back(0);	perm.push_back(1); perm.push_back(2);
		for(int i=1; i<permutations; ++i)
			std::next_permutation(perm.begin(), perm.end());
		//perm contains now the permuted indexes.
		for(int i=0; i<3; ++i)
		{
			for(int j=0; j<3; ++j)
				_axes[i][j] = axes[i][perm[j]];
			_axes[i] = _pDomain->_pLat[_map._mt]->toRotated(_axes[i]);
		}
		//normalize axes
		for(int i=0; i<3; ++i)
			_axes[i] /= _axes[i].abs();
		//if the specified vector is 0, then accept all orientations.
		if(edt->_doNotCreateIn.abs() == 0)
			break;
	}
	while(_axes[2]*edt->_doNotCreateIn == 0);
}


Cluster::~Cluster()
{
	while(_particles.size())
	{
		Particle *part = _particles.back();
		_pDomain->_pMesh->remove(part);
		delete part;
		_particles.pop_back();
	}
}

void Cluster::addTo(P_TYPE pt, ID &map, bool ivmodel)
{
	addMap(map, ClusterParam::pt2ID(map._mt, pt), ivmodel);
}

void Cluster::collapseMap(ID &map)
{  //careful! Changing the container while iterating on it.
	for(std::map<P_TYPE, unsigned>::iterator it=map._pt.begin(); it!=map._pt.end();)
		if(it->second == 0)
		{
			map._pt.erase(it->first);
			it = map._pt.begin();
		}
		else
			++it; //needed to avoid memory corruption and other problems.
	for(std::map<P_POS, unsigned>::iterator it=map._pos.begin(); it!=map._pos.end();)
		if(it->second == 0)
		{
			map._pos.erase(it->first);
			it = map._pos.begin();
		}
		else
			++it;
}

//takes a particle and adds it to the map and to the vector
//the particle is reproduced
Coordinates Cluster::addParticle(Kernel::SubDomain *pSub, Particle *p)
{
	P_TYPE fam= Domains::global()->PM()->getFamily(p->getPType());
	P_POS pos = Domains::global()->PM()->getPPos(p->getPType());
	P_TYPE mt1 = Domains::global()->PM()->getMaterial(_map._mt)._pt[0];
	P_TYPE mt2 = Domains::global()->PM()->getMaterial(_map._mt)._pt[1];
	bool binary = Domains::global()->PM()->getMaterial(_map._mt)._binary;
	vector<P_TYPE> fams; vector<P_POS> poss;
	if(pos == POS_V) //a pair, actually
	{
		poss.push_back(POS_0);
		poss.push_back((binary ? POS_1 : POS_0));
		fams.push_back(V_TYPE);
	}
	else
		poss.push_back(pos);
	fams.push_back(fam);
	bool ivmodel = _pDomain->_pClPar->getParams(_map._mt, _edtype)->_ivmodel;
	Kernel::Coordinates co = p->getCoordinates();
	while(fams.size())
	{
		fam = fams.back(); fams.pop_back();
		pos = poss.back(); poss.pop_back();
		//insert particle
		bool bIrreg = _pDomain->_pClPar->getParams(_map._mt, _edtype)->getAspect() == EDType::aspect_irregular;
		Kernel::MeshElement *pEle = _pElement;
		ID tempID = _map;
		addTo(p->getPType(), tempID, ivmodel);
		co = (bIrreg? p->getCoordinates() :
				_pDomain->_pClPar->getParams(_map._mt, _edtype)->surface(pSub->_rng, _axes, _center, tempID, p) + _center);
		if(_pDomain->_pMesh->checkMove(pEle, co) == Kernel::Mesh::JUMP_OK)
		{
			_particles.push_back(new Particle(fam, co, this, p->getOrig()));
			updateRadius(_particles.back(), +1);
			_pDomain->_pMesh->insert(_particles.back(), pEle);
		}
		else  //irregular whether we like or not
		{
			_particles.push_back(new Particle(fam, p->getCoordinates(), this, p->getOrig()));
			updateRadius(_particles.back(), +1);
			_pDomain->_pMesh->insert(_particles.back(), p->getElement());
		}
		if(!ivmodel)
		{
			//remove particles if needed
			P_POS annhPos = POS_V;
			if(fam == V_TYPE) annhPos = POS_I;
			if(fam == mt1) annhPos = POS_0;
			if(fam == mt2) annhPos = POS_1;
			if(annhPos != POS_V && _map._pos.find(annhPos) != _map._pos.end())
				removeParticle(fam);
			//
			P_TYPE annhPt = UNDEFINED_TYPE;
			if(pos == POS_I) annhPt = V_TYPE;
			if(pos == POS_0) annhPt = mt1;
			if(pos == POS_1) annhPt = mt2;
			if(annhPt != UNDEFINED_TYPE && _map._pt.find(annhPt) != _map._pt.end())
				removeParticle(annhPt);
		}
	}
	addTo(p->getPType(), _map, ivmodel); //updateMap
	return co;
}

//Everything in perfect state when calling this function!
//returns true if the cluster was killed by this function.
bool Cluster::fuzzModel(Kernel::SubDomain *pSub, P_TYPE pt, const Kernel::Coordinates &c, bool bKillIfInterface)
{
	P_TYPE fam= Domains::global()->PM()->getFamily(pt);
	const EDType *edPar = _pDomain->_pClPar->getParams(_map._mt, _edtype);
	const P_TYPE expand_impurity = edPar->_expand_impurity;
	bool bDilatation = fam == expand_impurity;
	bool touched_something = false;
	if(bKillIfInterface && expand_impurity != UNDEFINED_TYPE)
	//Sees if something was "touched"
	{
		float radius = edPar->_expand_capture_radius_P;
		if(edPar->_expand_capture_radius_M)
			radius += getRadius()*edPar->_expand_capture_radius_M;
		vector<Interface *> interfaces;
		_pDomain->_pMesh->fillInterfaces(c, radius, interfaces);
		for(vector<Interface *>::iterator it=interfaces.begin(); it!= interfaces.end(); ++it)
			if( Domains::global()->PM()->isGasLike((*it)->getElement(0)->getMaterial()) ||
				Domains::global()->PM()->isGasLike((*it)->getElement(1)->getMaterial()))
			{
				touched_something =  true;
				break;
			}
	}
	if(touched_something) //has to remove interface
	{
		float volume = edPar->_expand_impurity_volume_nm3;
		unsigned nImpurities = _map._pt.find(expand_impurity)->second;
		unsigned ratioAfter  = unsigned(volume* nImpurities   /_pElement->getVolume());
		Kernel::Mesh * pMesh = _pDomain->_pMesh;
		Kernel::MeshElement *pEle = _pElement;
		_pDomain->_pRM->remove(this, _pElement);
		delete this;
		while(ratioAfter--)
			pMesh->changeNumberOfCells(pSub, pEle, -1);
		return true;
	}
	if(bDilatation)
	{
		float volume = _pDomain->_pClPar->getParams(_map._mt, _edtype)->_expand_impurity_volume_nm3;
		unsigned nImpurities = _map._pt.find(expand_impurity)->second;
		unsigned ratioAfter  = unsigned(volume* nImpurities   /_pElement->getVolume());
		unsigned ratioBefore = unsigned(volume*(nImpurities-1)/_pElement->getVolume());
		if(ratioAfter == ratioBefore+1)
			_pDomain->_pMesh->changeNumberOfCells(pSub, _pElement, 1);//expand the material by 1 volume element.
	}
	return false;
}

// performing actions so far are
// to emit the particle.
// migration
void Cluster::perform (Kernel::SubDomain *pSub, unsigned eventType)
{
	pSub->_evLog.performed(_map._mt, Event::CLUSTER, _edtype, _hash, eventType, _state);
	if(eventType < MAX_PARTICLES)
		emit(pSub, eventType);
	else if(eventType == MAX_PARTICLES)
		migrate(pSub);
	else if(eventType == MAX_PARTICLES +1)
		transform(pSub, _pDomain->_pClPar->getParams(_map._mt, _edtype)->_transformTo);
	else if(eventType == MAX_PARTICLES+2)
		transform(pSub, _pDomain->_pClPar->getParams(_map._mt, _edtype)->_transformFrom);
	else
		recombine(pSub);
}

/*
 * It removes the last particle incorporated to try to keep the shape OK
 * But the particle is emitted from a position that depends on the shape.
 */
void Cluster::emit(Kernel::SubDomain *pSub, unsigned eventType)
{
	_pDomain->_pRM->remove(this,_pElement);
	const EDType *edt = _pDomain->_pClPar->getParams(_map._mt, _edtype);
	ID oldMap = _map;
	emitFrom(eventType, _map); // Update my map.
	Particle const *toEmit = edt->emitFrom(pSub->_rng, _center, _particles);
	Kernel::MeshElement *pEle = toEmit->getElement();
	Coordinates toEmitCoord = toEmit->getCoordinates(), temp_c = toEmit->getCoordinates();
	const EDType::CLType *otherCl = _pDomain->_pClPar->getParams(_edtype, _map);
	P_TYPE otherPt = ClusterParam::ID2pt(_map);
	if(otherPt != UNDEFINED_TYPE)
		breakPosition(pSub, temp_c, pEle, _pDomain->_pMPPar->_interact_radius[_map._mt][otherPt][eventType]);
	else
		breakPosition(pSub, temp_c, pEle, otherCl->_interactMP_radius[eventType]);
	new MobileParticle(pSub, eventType, 0,  _pDomain, temp_c, toEmit->getOrig()); //emit one particle

	P_TYPE fam = Domains::global()->PM()->getFamily(eventType);
	P_POS  pos = Domains::global()->PM()->getPPos(eventType);
	P_TYPE mt1 = Domains::global()->PM()->getMaterial(_map._mt)._pt[0];
	P_TYPE mt2 = Domains::global()->PM()->getMaterial(_map._mt)._pt[1];
	bool binary = Domains::global()->PM()->getMaterial(_map._mt)._binary;
	vector<P_TYPE> fams; vector<P_POS> poss;
	if(pos == POS_V) //a pair, actually
	{
		poss.push_back(POS_0);
		poss.push_back((binary ? POS_1 : POS_0));
		fams.push_back(V_TYPE);
	}
	else
		poss.push_back(pos);
	fams.push_back(fam);
	while(fams.size())
	{
		fam = fams.back(); fams.pop_back();
		pos = poss.back(); poss.pop_back();
		//pos
		if(oldMap._pos.find(pos) == oldMap._pos.end())  //increase it as family
		{
			P_TYPE newPt = UNDEFINED_TYPE;
			if(pos == POS_I) newPt = V_TYPE;
			if(pos == POS_0) newPt = mt1;
			if(pos == POS_1) newPt = mt2;
			assert(newPt != UNDEFINED_TYPE);
			Particle *pPart = new Particle(newPt, toEmitCoord, this, toEmitCoord);
			_pDomain->_pMesh->insert(pPart, 0);
			_particles.push_back(pPart);
			updateRadius(pPart, +1);
		}
		if(oldMap._pt.find(fam) != oldMap._pt.end())
			removeParticle(fam);
	}
	vector<Particle *> dummy;
	checkSize(pSub, true, dummy);
}

bool Cluster::emitFrom(P_TYPE pt, ID &map)
{
	P_TYPE fam = Domains::global()->PM()->getFamily(pt);
	P_POS  pos = Domains::global()->PM()->getPPos(pt);
	P_TYPE mt1 = Domains::global()->PM()->getMaterial(map._mt)._pt[0];
	P_TYPE mt2 = Domains::global()->PM()->getMaterial(map._mt)._pt[1];
	bool binary = Domains::global()->PM()->getMaterial(map._mt)._binary;
	vector<P_TYPE> fams; vector<P_POS> poss;
	if(pos == POS_V) //a pair, actually
	{
		poss.push_back(POS_0);
		poss.push_back((binary? POS_1:POS_0));
		fams.push_back(V_TYPE);
	}
	else
		poss.push_back(pos);
	fams.push_back(fam);
	while(fams.size())
	{
		fam = fams.back(); fams.pop_back();
		pos = poss.back(); poss.pop_back();
		//fam
		if(map._pt.find(fam) != map._pt.end())
			map._pt[fam]--;
		else //increase it as position
		{
			if(fam != mt1 && fam != mt2 && fam != V_TYPE)
				return false;
			P_POS newPos = (fam == V_TYPE? POS_I : (fam == mt1? POS_0 : POS_1));
			if(map._pos.find(newPos) == map._pos.end())
				map._pos[newPos] = 0;
			map._pos[newPos]++;
		}
		//pos
		if(map._pos.find(pos) != map._pos.end())
			map._pos[pos]--;
		else //increase it as family
		{
			P_TYPE newPt = UNDEFINED_TYPE;
			if(pos == POS_I) newPt = V_TYPE;
			if(pos == POS_0) newPt = mt1;
			if(pos == POS_1) newPt = mt2;
			assert(newPt != UNDEFINED_TYPE);
			if(map._pt.find(newPt) == map._pt.end())
				map._pt[newPt] = 0;
			map._pt[newPt]++;
		}
	}
	collapseMap(map);
	return true;
}

void Cluster::migrate(Kernel::SubDomain * pSub)
{
	assert(_particles.size());
	Kernel::MeshElement *pEle = _pElement;
	const Coordinates orig = _center;
	int random = int(12.* pSub->_rng.rand());
	int ax, sign = (random%2)*2-1;
	const EDType *edType = _pDomain->_pClPar->getParams(_map._mt, _edtype);
	switch(edType->_migration)
	{
	case EDType::mig_perpendicular:
		ax = 2;
		break;
	case EDType::mig_parallel:
		ax = random/6;
		break;
	case EDType::mig_3d:
		ax = random/4;
		break;
	default:
		ax = 0;
		break;
	}
	Coordinates delta = _axes[ax]*sign* edType->_lambda, dummy(0,0,0);
	Kernel::Mesh::JUMP_ACTIONS jmp = _pDomain->_pMesh->jumpPosition(_center, delta, pEle, 0, dummy, 1);
	if (jmp != Kernel::Mesh::JUMP_OK)
	{
		_center = orig;
		return;
	}
	if(pEle->getMaterial() != _map._mt) //interface
	{
		_center = orig;
		if(interactSurface(pSub, _pElement, pEle) == NO_INTERFACE)
			ERRORMSG("Cluster: Could not find interface to react to");
		return;
	}

	if (_pElement->getBAtoms() != pEle->getBAtoms()) //change in formation and migrations
	{
		float kT = _pDomain->_pRM->getkT();
		// Segregation and concentrations
		IO::Arrhenius mig = _myParams->_arr[0](_pElement);
		IO::Arrhenius form =_myParams->_eForm(_pElement);
		IO::Arrhenius newmig = _myParams->_arr[0](pEle);
		IO::Arrhenius newform = _myParams->_eForm(pEle);

		double diff = mig._ener + form._ener;
		double newdiff = newmig._ener + newform._ener;
		if ( newdiff > diff)
		{
			double prob = (newform._pref / form._pref) * (newmig._pref / mig._pref) * std::exp( - (newdiff - diff) / kT);
			if (pSub->_rng.rand() > prob)
			{
				_center = orig;
				return;
			}
		}
	}
	Kernel::Mesh::ClusterInfo tmp;
	vector<Kernel::Mesh::ClusterInfo> after;
	for(vector<Particle *>::iterator it=_particles.begin(); it!=_particles.end(); ++it)
	{
		tmp._coord = (*it)->getCoordinates();
		tmp._orig  = (*it)->getOrig();
		tmp._pEle  = (*it)->getElement();
		tmp._pPart = (*it);
		if(_pDomain->_pMesh->jumpPosition(tmp._coord, delta, tmp._pEle, 0, tmp._orig, 1) != Kernel::Mesh::JUMP_OK)
		{
			_center = orig;
			return;
		}
		if(tmp._pEle->getMaterial() != _map._mt) //interface
		{
			if(interactSurface(pSub, (*it)->getElement(), tmp._pEle) == NO_INTERFACE)
				ERRORMSG("Cluster: Could not find interface to react to");
			_center = orig;
			return;
		}
		after.push_back(tmp);
	}
	if(_pElement != pEle) //careful here, this function can be slow...
	{
		if(_pDomain->_pRM->getSubDomains() != 1 || _pElement->getBAtoms() != pEle->getBAtoms())
		{
			_pDomain->_pRM->remove(this, _pElement); //to udpate subdomains and levels.
			_pElement = pEle;
			_pDomain->_pRM->insert(this, pEle);
		}
		else
			_pElement = pEle;
	}
	for(vector<Kernel::Mesh::ClusterInfo>::iterator it=after.begin(); it!=after.end(); ++it)
	{
		if(it->_pEle != it->_pPart->getElement())
		{
			_pDomain->_pMesh->remove(it->_pPart);
			it->_pPart->setCoordinates(it->_coord);
			it->_pPart->setOrig(it->_orig);
			_pDomain->_pMesh->insert(it->_pPart, it->_pEle);
		}
		else
			it->_pPart->setCoordinates(it->_coord);
	}
	//All particles out, but updated
	//Look for neighbors: Caution! This function can be VERY slow if many particles are around
	vector<Particle *> interact;
	bool lookFor = true;
	if(_particles.size() > 4)
	{
		vector<Kernel::MeshElement *> elementsToLook;
		const float capture = Domains::global()->PM()->getMaximumCaptureRadius(_map._mt);
		_pDomain->_pMesh->fillNeighborsOneMat(_pElement, _center, _radius + capture, elementsToLook);
		unsigned totalParticles = 0;
		for(vector<Kernel::MeshElement *>::iterator it=elementsToLook.begin(); it!=elementsToLook.end(); ++it)
			totalParticles += (*it)->getNumberOfParticles();
		if(totalParticles == _particles.size())
			lookFor = false;
		if(totalParticles < _particles.size())
			ERRORMSG(capture << " " << _radius << " " << getRadius() <<	" " << _center << " " << pEle->getCoords(0.5, 0.5, 0.5) << std::endl <<
					" " << elementsToLook.size() << " elements " << std::endl <<
			        "Cluster::migrate " << totalParticles << " cannot be < than " << _particles.size());
	}
	if(lookFor)
		for(vector<Kernel::Mesh::ClusterInfo>::iterator it=after.begin(); it!=after.end(); ++it)
			_pDomain->_pMesh->getInteractions(pSub, it->_pPart, interact);
	if(interact.size())
	{
		Particle *intPar = interact[int(pSub->_rng.rand() * interact.size())];
		std::vector<Particle *> dummy;
		intPar->getDefect()->interact(pSub, this, intPar, dummy);
	}
}

void Cluster::transform(Kernel::SubDomain *pSub, unsigned edT)
{
	assert(edT != _edtype);
	_pDomain->_pRM->remove(this, _pElement);
	Cluster *pMC = new Cluster(pSub, _pDomain, edT, _center);
	pMC = pMC->growDefect(pSub, _map);
	if(pMC)
		_pDomain->_pRM->insert(pMC, pMC->_pElement);
	delete this;
}

void Cluster::recombine(Kernel::SubDomain *pSub)
{
	bool bRemoved = false;
	_pDomain->_pRM->remove(this, _pElement);
	map<P_POS,  unsigned>::iterator itPos = _map._pos.find(POS_I);
	map<P_TYPE, unsigned>::iterator itPt  = _map._pt.find(V_TYPE);
	if(itPos != _map._pos.end() && itPt != _map._pt.end())
	{  //remove an IV pair
		itPos->second--;
		itPt->second--;
		removeParticle(V_TYPE);
		bRemoved = true;
	}
	P_POS pos[2] = { POS_0, POS_1 };
	for(unsigned p=0; p<2; ++p)
	{
		itPos = _map._pos.find(pos[p]);
		P_TYPE pt = Domains::global()->PM()->getMaterial(_map._mt)._pt[p];
		itPt  = _map._pt.find(pt);
		if(itPos != _map._pos.end() && itPt != _map._pt.end())
		{
			itPos->second--;
			itPt->second--;
			removeParticle(pt);
			if(bRemoved) //IV or SiC + CSi
				break;
			bRemoved = true;

		}
	}
	if(!bRemoved)
		WARNINGMSG("Cluster " << name() << " has a non-zero recombination rate, but nothing to recombine!");
	collapseMap(_map);
	vector<Particle *> dummy;
	checkSize(pSub, true, dummy);
}

bool Cluster::recombination(ID &theMap)
{
	bool bRemoved = false;
	map<P_POS,  unsigned>::iterator itPos = theMap._pos.find(POS_I);
	map<P_TYPE, unsigned>::iterator itPt  = theMap._pt.find(V_TYPE);
	if(itPos != theMap._pos.end() && itPt != theMap._pt.end())
	{  //remove an IV pair
		itPos->second--;
		itPt->second--;
		bRemoved = true;
	}
	P_POS pos[2] = { POS_0, POS_1 };
	for(unsigned p=0; p<2; ++p)
	{
		itPos = theMap._pos.find(pos[p]);
		P_TYPE pt = Domains::global()->PM()->getMaterial(theMap._mt)._pt[p];
		itPt  = theMap._pt.find(pt);
		if(itPos != theMap._pos.end() && itPt != theMap._pt.end())
		{
			itPos->second--;
			itPt->second--;
			if(bRemoved) //IV pair or another recombination.
				break;
			bRemoved = true;
		}
	}
	collapseMap(theMap);
	return bRemoved;
}

void Cluster::removeParticle(P_TYPE fam)
{
	vector<Particle *>::reverse_iterator rit=_particles.rbegin();
	for(; rit!=_particles.rend(); ++rit)
		if((*rit)->getPType() == fam) //it found the particle to annihilate
		{
			_pDomain->_pMesh->remove(*rit);
			delete *rit;
			_particles.erase(--(rit.base()));
			updateRadius(*rit, -1);
			break;
		}
	assert(rit != _particles.rend());
}

float Cluster::getRate(unsigned eventType, float kT) const
{
	if(eventType < MAX_PARTICLES)
		return emissionRate(getElement(), _edtype, _myParams, _map, eventType, kT);
	else  if(eventType == MAX_PARTICLES + 3)  //recombination
		return recombinationRate(getElement(), _edtype, _myParams, _map,  kT);
	else //Migration, transformations.
		return _myParams->_arr[eventType - MAX_PARTICLES](getElement()).getRate(kT);
}

double Cluster::recombinationRate(const Kernel::MeshElement *pEle, unsigned edtype, const EDType::CLType *myParams, const ID &thisMap, double kT)
{
	const Kernel::Domain *pDomain = pEle->getDomain();
	double rate = myParams->_arr[3](pEle).getRate(kT);
	if(rate == 0)
		return 0;
	ID theMap(thisMap);
	//TODO: if the pos_0 cannot be recombined (is not defined) but pos_1 can, this does not work properly.
	if(!recombination(theMap))
		return 0;
	float diffForm = 0;
	const EDType::CLType *otherCl = pDomain->_pClPar->getParams(edtype, theMap);
	P_TYPE pt = ClusterParam::ID2pt(theMap);
	if(theMap._pt.size() == 0 || pt != UNDEFINED_TYPE) //the remaining cluster when it is a particle
	{
		if(pt != UNDEFINED_TYPE)
			diffForm += pDomain->_pMPPar->_form[thisMap._mt][pt](pEle)._ener;
	}
	else
	{
		if(otherCl == 0)
		{
			WARNINGMSG("Cluster " << Domains::global()->PM()->getIDName(thisMap) << " recombination rate is non null but cluster does not exist");
			return 0;
		}
		else
			diffForm += otherCl->_eForm(pEle)._ener;
	}
	//IV recombinations or emissions included in formation!
	diffForm -= myParams->_eForm(pEle)._ener;
	float actEner = std::max(diffForm, float(0.));
	//diffpot plus migration
	return rate*exp(-actEner / kT);
}

double Cluster::emissionRate(const Kernel::MeshElement *pEle, unsigned edtype, const EDType::CLType *myParams, const ID &thisMap, unsigned eventType, double kT)
{
	if(eventType >= Domains::global()->PM()->getNParticles() ||
			eventType < Domains::global()->PM()->getNFamilies()) //non defined particles
		return 0;
	Kernel::Domain *pDomain = pEle->getDomain();
	Kernel::M_TYPE mt = thisMap._mt;
	if(myParams->_pref[eventType] == 0)
		return 0;
	ID theMap(thisMap);

	bool bPossible = emitFrom(eventType, theMap);
	if(!bPossible)
	{
		WARNINGMSG("Cluster " <<
				Domains::global()->PM()->getIDName(thisMap) << " prefactor for particle " <<
				Domains::global()->PM()->getParticleName(mt, eventType) << " is non null (" <<
				myParams->_pref[eventType] << "), but cluster cannot emit such particle");
		return 0;
	}
	const EDType::CLType *otherCl = pDomain->_pClPar->getParams(edtype, theMap);
	P_TYPE pt = ClusterParam::ID2pt(theMap);
	float diffForm = pDomain->_pMPPar->_form[mt][eventType](pEle)._ener;
	IO::Arrhenius prodMig;
	if(pt != UNDEFINED_TYPE) //the remaining cluster when it is a particle
	{
		if(pDomain->_pMPPar->_canInteract[mt][eventType][pt] == false)
			return 0; //microsc. reversibility.
		//no prodMig, to avoid double emission when calling with the rate of the other emitting particle
		diffForm += pDomain->_pMPPar->_form[mt][pt](pEle)._ener;
	}
	else
	{
		if(otherCl == 0)
		{
			WARNINGMSG("Cluster " <<
					Domains::global()->PM()->getIDName(thisMap) << " prefactor for particle " <<
					Domains::global()->PM()->getParticleName(mt, eventType) << " is non null (" <<
					myParams->_pref[eventType] << "), but cluster does not exist");
			return 0;
		}
		else
		{
			if(otherCl->_interactMP[eventType] == false) //microscopical reversibility
				return 0;
			diffForm += otherCl->_eForm(pEle)._ener;
			prodMig   = otherCl->_arr[0](pEle);
		}
	}
	//IV recombinations or emissions included in formation!
	diffForm -= myParams->_eForm(pEle)._ener;
	//1. particle moves; 2. Cluster moves
	float actEner  = std::max(diffForm, float(0.)) + pDomain->_pMPPar->_arr[mt][eventType][0][0](pEle)._ener;
	float actEner2 = std::max(diffForm, float(0.)) + prodMig._ener;
	//diffpot plus migration
	int f1 = (pDomain->_pMPPar->_arr[mt][eventType][0][0](pEle)._pref != 0 ? 1: 0);
	int f2 = prodMig._pref != 0? 1: 0;
	return myParams->_pref[eventType] * ( f1*exp(-actEner / kT) + f2*exp(-actEner2 / kT) );
}

//the barrier is considered in canInteract
// the vector of particles is not implemented yet, because it is used in Interface only.
Defect * Cluster::interact(Kernel::SubDomain *pSub, Defect *def, Particle *, std::vector<Particle *> &parts)
{
	_pDomain->_pRM->remove(this, _pElement);

	if(def->getEType() == MOBILEPARTICLE) //MC + MP
	{
		MobileParticle *mp = static_cast<MobileParticle *>(def);
		pSub->_reLog.reaction(_map._mt, CLUSTER, getEDType(), _state, MOBILEPARTICLE, mp->getPType(),mp->getState());
		if(_myParams->_interactMP_sink[mp->getPType()])
		{
			delete this;
			return 0;
		}
		Kernel::Coordinates co = addParticle(pSub, mp);
		P_TYPE thePt = mp->getPType();
		delete mp;
		Defect *def = checkSize(pSub, true, parts);
		if(def && fuzzModel(pSub, thePt, co, true))
			return 0;
		if(def && _pDomain->_pClPar->getParams(_map._mt, _edtype)->_percolation) //Cluster + cluster
		{
			vector<Particle *> interact;
			if(_particles.size())
				_pDomain->_pMesh->getInteractions(pSub, _particles.back(), interact);
			if(interact.size())
			{
				Particle *intPar = interact[int(pSub->_rng.rand() * interact.size())];
				std::vector<Particle *> dummy;
				return intPar->getDefect()->interact(pSub, this, intPar, dummy);
			}
		}
		return def;
	}
	if(def->getEType() == CLUSTER) //MC + MC
    {
		Cluster *mc = static_cast<Cluster *>(def);
		assert(def != this);
		_pDomain->_pRM->remove(mc, mc->_pElement);
		pSub->_reLog.reaction(_map._mt, CLUSTER, getEDType(), _state, CLUSTER, mc->getEDType(), mc->getState());
		unsigned deft0 = _edtype, deft1 = mc->_edtype;
		unsigned sizet0 = _particles.size(), sizet1 = mc->_particles.size();
		int edType_temp = _pDomain->_pClRePar->_interactions[_map._mt][deft0][deft1](sizet0, sizet1);
		if (edType_temp == -1)
			ERRORMSG("Error in cluster interactions for types " << deft0 << " " << deft1 << " sizes " << sizet0 << " " << sizet1);
		unsigned edType = unsigned(edType_temp);
		vector<std::pair<Kernel::MeshElement *, P_TYPE> > fuzzElements; //temporal storage to avoid calling fuzz at the middle of a
			// non well defined state...

		Coordinates newCenter;
		_pDomain->_pMesh->middle(newCenter, _center, mc->_center);
		Cluster *mcNew = new Cluster(pSub, _pDomain, edType, newCenter);
		bool ivmodel = _pDomain->_pClPar->getParams(_map._mt, edType)->_ivmodel;
		ID newMap = _map;
		addMap(newMap, mc->_map, ivmodel);
		delete mc;
		delete this;
		// grow calls fuzz, and fuzz requires everything in perfect state!
		mcNew = mcNew->growDefect(pSub, newMap);
		if(mcNew)
			mcNew->_pDomain->_pRM->insert(mcNew, mcNew->_pElement);
		return mcNew; //no check for mcNew, grow defect checks it.

	}
	ERRORMSG("Cluster: Interaction not allowed!");
	return 0;
}

//adds ID2 into ID1
void Cluster::addMap(ID &ID1, const ID &ID2, bool ivmodel)
{
	P_TYPE mt1 = Domains::global()->PM()->getMaterial(ID1._mt)._pt[0];
	P_TYPE mt2 = Domains::global()->PM()->getMaterial(ID1._mt)._pt[1];

	for(map<P_TYPE, unsigned>::const_iterator it=ID2._pt.begin(); it != ID2._pt.end(); ++it)
		for(unsigned howMany=0; howMany < it->second; ++howMany)
		{
			P_POS annhPos = POS_V;
			if(it->first == V_TYPE) annhPos = POS_I;
			if(it->first == mt1)    annhPos = POS_0;
			if(it->first == mt2)    annhPos = POS_1;
			if(!ivmodel && annhPos != POS_V && ID1._pos.find(annhPos) != ID1._pos.end())
			{
				if(--ID1._pos[annhPos] == 0)
					ID1._pos.erase(annhPos);
			}
			else
			{
				if(ID1._pt.find(it->first) == ID1._pt.end())
					ID1._pt[it->first] = 0;
				ID1._pt[it->first]++;
			}
		}
	for(map<P_POS, unsigned>::const_iterator it=ID2._pos.begin(); it != ID2._pos.end(); ++it)
		for(unsigned howMany=0; howMany < it->second; ++howMany)
		{
			P_TYPE annhPt = UNDEFINED_TYPE;
			if(it->first == POS_I) annhPt = V_TYPE;
			if(it->first == POS_0) annhPt = mt1;
			if(it->first == POS_1) annhPt = mt2;
			if(!ivmodel && annhPt != UNDEFINED_TYPE && ID1._pt.find(annhPt) != ID1._pt.end())
			{
				if(--ID1._pt[annhPt] == 0)
					ID1._pt.erase(annhPt);
			}
			else
			{
				if(ID1._pos.find(it->first) == ID1._pos.end())
					ID1._pos[it->first] = 0;
				ID1._pos[it->first]++;
			}
		}
	collapseMap(ID1);
}

//to be called with the defect removed from Scheduler
Defect * Cluster::checkSize(Kernel::SubDomain *pSub, bool bInsert, 	vector<Particle *> &parts)
{
	if(_particles.size() == 0)
	{
		delete this;
		return 0;
	}
	if(ClusterParam::ID2pt(_map) != UNDEFINED_TYPE)
	{
		Particle *otherEmit = _particles.back();
		for(vector<Particle *>::iterator it=_particles.begin(); it!=_particles.end(); ++it)
			if(Domains::global()->PM()->getIorV(_map._mt, (*it)->getPType()) == Kernel::IV_NONE)
				otherEmit = *it;
		MobileParticle *mp = new MobileParticle(pSub, ClusterParam::ID2pt(_map),0,
				_pDomain, otherEmit->getCoordinates(), otherEmit->getOrig());
		delete this;
		return 0;
	}
	else
	{
		if(_particles.front()->getElement() != _pElement)
		{
			_pElement = _particles.front()->getElement();
			_center   = _particles.front()->getCoordinates();
		}
		_myParams = getMyParam();
		if(bInsert)
			_pDomain->_pRM->insert(this, _pElement);
		return this;
	}
}

const EDType::CLType * Cluster::getMyParam()
{
	Kernel::M_TYPE mt = _particles.front()->getElement()->getMaterial();
	_hash = _pDomain->_pClPar->getParams(mt, _edtype)->_hash.cluster2hash(_map);
	const EDType::CLType *data = _pDomain->_pClPar->getParams(mt, _edtype, _hash);
	if(!data)
		ERRORMSG("Cluster::getMyParam: " << Domains::global()->PM()->getIDName(_map) << " is not defined");
	return data;
}

bool Cluster::canInteract(Kernel::SubDomain *pSub, const Particle *dyn, const Particle *sta) const
{
	Defect *def = dyn->getDefect();
	float eForm_init = 0;
	float eForm_end = 0;
	if(def->getState()!=0 || _state !=0)
		return false;
	if(Domains::global()->PM()->isAmorphous(_pElement->getMaterial()) || Domains::global()->PM()->isAmorphous(def->getElement()->getMaterial()))
		return false;
	float distSq = 0;
	if(def->getEType() == MOBILEPARTICLE)
	{
		P_TYPE pt = dyn->getPType();
		if(!_myParams->_interactMP[pt])
			return false;
		distSq = _myParams->_interactMP_radius[pt] * _myParams->_interactMP_radius[pt];
		{	//check distances
			Coordinates n = dyn->getCoordinates();
			_pDomain->_pMesh->setPeriodicRelative(sta->getCoordinates(), n);
			if(n*n >= distSq)
				return false;
		}
		if(_myParams->_interactMP_sink[pt])
		{
			eForm_end = eForm_init = 0; //sink interaction always happens.
		}
		else
		{
			ID theMap(_map);
			bool ivmodel = _pDomain->_pClPar->getParams(_map._mt, _edtype)->_ivmodel;
			addTo(pt, theMap, ivmodel);
			const EDType::CLType *otherCl = _pDomain->_pClPar->getParams(_edtype, theMap);
			//now the mobile particle
			float eForm_binding = _pDomain->_pMPPar->_form[_map._mt][pt](_pElement)._ener;
			P_TYPE ptEnd = ClusterParam::ID2pt(theMap);
			if(ptEnd != UNDEFINED_TYPE)
				eForm_end = _pDomain->_pMPPar->_form[_map._mt][ptEnd](_pElement)._ener;
			else
			{
				if(otherCl == 0)
				{
					WARNINGMSG("Cluster::canInteract produces a cluster " <<
							Domains::global()->PM()->getIDName(theMap) << "that is NOT defined!");
					return false;
				}
				eForm_end = otherCl->_eForm(_pElement)._ener;
			}
			eForm_init += _myParams->_eForm (_pElement)._ener + eForm_binding;
		}
	}
	else if(def->getEType() == CLUSTER)
	{
		const Cluster *pMC = static_cast<const Cluster *>(def);
		const EDType *data = _pDomain->_pClPar->getParams(_map._mt, _edtype);
		std::map<unsigned, float>::const_iterator it = data->_interactMC.find(pMC->_edtype);
		if(it == data->_interactMC.end())
			return false;
		distSq = it->second * it->second;
		{	//check distances
			Coordinates n = dyn->getCoordinates();
			_pDomain->_pMesh->setPeriodicRelative(sta->getCoordinates(), n);
			if(n*n >= distSq)
				return false;
		}
		unsigned deft1 = _edtype, deft0 = pMC->_edtype;
		unsigned sizet1 = _particles.size(), sizet0 = pMC->_particles.size();
		int edType_temp = _pDomain->_pClRePar->_interactions[_map._mt][deft0][deft1](sizet0, sizet1);
		if (edType_temp == -1)
			ERRORMSG("Error in cluster interactions for types " << deft0 << " " << deft1 << " sizes " << sizet0 << " " << sizet1);
		unsigned edType = unsigned(edType_temp);
		ID newMap = _map;
		bool ivmodel = _pDomain->_pClPar->getParams(_map._mt, edType)->_ivmodel;
		addMap(newMap, pMC->_map, ivmodel);
		P_TYPE ppt = ClusterParam::ID2pt(newMap);
		if(newMap._pt.size())
		{
			if(ppt == UNDEFINED_TYPE)
			{
				const EDType::CLType * cltype = _pDomain->_pClPar->getParams(edType, newMap);
				if(!cltype)
					return false; //cluster does not exist, it cannot react top form it. (i.e. I2 + I9 with a maximum I10)
				eForm_end = cltype->_eForm(_pElement)._ener;
			}
			else
				eForm_end = _pDomain->_pMPPar->_form[_map._mt][ppt](_pElement)._ener;
		}
		//compute probs.
		eForm_init += pMC->_myParams->_eForm(_pElement)._ener + _myParams->_eForm(_pElement)._ener;
	}
	else if(def->getEType() == INTERFACE)
		return def->canInteract(pSub, sta, dyn);
	else
		ERRORMSG("Cluster. Reaction not allowed with defect type " << def->getEType());
	double arg = eForm_end - eForm_init;
	if(def->getEType() == MOBILEPARTICLE)  //metastable
		return (arg < 0? true : pSub->_rng.rand() < exp(arg/_pDomain->_pRM->getkT()));
	else
		return (arg <= 0);
}

double Cluster::getRadius() const
{
	//compute radius
	const Coordinates first = _particles[0]->getCoordinates();
	Coordinates cdm(0,0,0);
	double myradius = 0;
	unsigned howMany = 1;
	for(vector<OKMC::Particle *>::const_iterator itPart = _particles.begin()+1; itPart!=_particles.end(); ++itPart)
	{
		Coordinates second = (*itPart)->getCoordinates();
		_pDomain->_pMesh->setPeriodicRelative(first, second);
		cdm += second;
		howMany++;
	}
	cdm /= howMany;
	cdm += first;
	double max2 = 0;

	for(vector<OKMC::Particle *>::const_iterator itPart = _particles.begin()+1; itPart!=_particles.end(); ++itPart)
	{
		Coordinates second = (*itPart)->getCoordinates();
		_pDomain->_pMesh->setPeriodicRelative(cdm, second);
		if(second*second > max2)
		{
			max2 = second*second;
			myradius = std::sqrt(max2);
		}
	}
	return myradius;
}

//the whole cluster is removed, not the best thing, but the easiest!
void Cluster::deletePart(Kernel::SubDomain *pSub, std::vector<Particle *> &parts)
{
	_pDomain->_pRM->remove(this,_pElement);
	//dissolve the whole thing
	for(vector<Particle *>::iterator it=_particles.begin(); it!=_particles.end(); ++it)
	{
		if(*it != parts.back())
		{
			Kernel::M_TYPE mtCheck = (*it)->getElement()->getMaterial();
			for(std::map<P_POS, unsigned>::iterator itp = _map._pos.begin(); itp != _map._pos.end(); ++itp)
			{
				P_TYPE ptCheck = Domains::global()->PM()->getParticle(mtCheck, (*it)->getPType(), itp->first);
				if(ptCheck != UNDEFINED_TYPE && Domains::global()->PM()->isParticleDefined(ptCheck,mtCheck))
				{
					itp->second--;
					if(itp->second == 0)
						_map._pos.erase(itp->first);
					MobileParticle *mp = new MobileParticle(pSub, ptCheck, 0, _pDomain, (*it)->getCoordinates(), (*it)->getOrig());
					std::replace(parts.begin(), parts.end(), *it, static_cast<Particle *>(mp));
					break;
				}
			}
		}
		if(std::remove(parts.begin(), parts.end(), *it) != parts.end())
			parts.pop_back();
		_pDomain->_pMesh->remove(*it);
		delete *it;
	}
	_particles.clear();
	delete this;
}

void Cluster::updateRadius(const Particle *part, int adding)
{
	if(adding > 0)
	{
		Coordinates c = part->getCoordinates();
		_pDomain->_pMesh->setPeriodicRelative(_center, c);
		double myRad = sqrt(c*c);
		if(myRad > _radius)
			_radius = myRad;
		return;
	}
	//substrating particle
	_radius = 0;
	for(vector<Particle *>::iterator it=_particles.begin(); it!=_particles.end(); ++it)
		updateRadius(*it, +1);
}

}

