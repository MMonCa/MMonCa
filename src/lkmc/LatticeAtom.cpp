/**
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

#include "LatticeAtom.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "kernel/RateManager.h"
#include "domains/MCClient.h"
#include "io/Diagnostic.h"
#include "io/SaveCmd.h"
#include "LatticeParam.h"
#include "LKMCModel.h"
#include <cstdlib>

using Kernel::Mesh;
using Kernel::Coordinates;
using std::vector;

namespace LKMC {

int LatticeAtom::_num =0;

LatticeAtom::LatticeAtom(LASTATE state, Kernel::Domain *p, Kernel::M_TYPE basicMt, Kernel::P_TYPE type, const Coordinates &c) :
	Event(p), LatticeSite(basicMt, type, c), _state(state), _maxNeigh(0)

{
	_number = _num++;
	_neighbors = 0;

	if (_pDomain->_pLaPar[basicMt]->_mapToGrid)
		_pDomain->_pMesh->getParticleToNodeHandler()->insert(this);
}

LatticeAtom::LatticeAtom(std::istream &is) : Event(is), LatticeSite(is), _maxNeigh(0)
{
	_number = _num++;
	int state;
	is >> state;
	_state = LASTATE(state);
	_neighbors = 0;

	if (_pDomain->_pLaPar[_basicMat]->_mapToGrid)
		_pDomain->_pMesh->getParticleToNodeHandler()->insert(this);
}

void LatticeAtom::restart(std::ostream &os) const
{
	Event::restart(os);
	LatticeSite::restart(os);
	os << _state << " ";
}

LatticeAtom::~LatticeAtom()
{
	_pDomain->_pMesh->getParticleToNodeHandler()->remove(this);

	for(unsigned i = 0; i < _maxNeigh; ++i)
		if(_neighbors[i])
			_neighbors[i]->removeNeighbor(this);
	delete [] _neighbors;
}

void LatticeAtom::setState(LASTATE st) 
{
	_pElement->decLA(getPerformed()); 
	_state = st; 
	_pElement->incLA(getPerformed()); 
}

void LatticeAtom::error(LatticeAtom * other, char const * const name) {
	std::multimap<std::string, Coordinates> result;
	LOWMSG(_number << " Orig atom is " << _coord);
	result.insert(std::pair<std::string, Coordinates>("original", _coord));
	LOWMSG(other->_number << " inserted atom is " << other->_coord);
	result.insert(std::pair<std::string, Coordinates>("inserted", other->_coord));
	std::string level("neighbor1");
	for(unsigned i=0; i < _maxNeigh; ++i)
	{
		if(i == first()) {
			level = "neighbor2";
		}
		else if(i == second()) {
			level = "neighbor3";
		}
		if(_neighbors[i])
		{
			Kernel::Coordinates ci = _neighbors[i]->getCoordinates();
			Kernel::Coordinates cj = getCoordinates();
			_pDomain->_pMesh->setPeriodicRelative(ci, cj);
			result.insert(std::pair<std::string, Coordinates>(level, _neighbors[i]->getCoordinates()));
			LOWMSG(_neighbors[i]->_number << " Neighbor " << i << " "<< _neighbors[i]->getCoordinates() <<
					" " << sqrt(cj._x*cj._x + cj._y*cj._y +	cj._z*cj._z) );
		}
	}
	std::string filename("neig-error-");
	std::ofstream out(filename + name + ".xyz");
	out << result.size() << "\n\n";
	for(auto const& item : result) {
		out << item.first << ' ' << item.second._x << ' ' << item.second._y << ' ' << item.second._z << '\n';
	}
}

void LatticeAtom::insertNeighbors()
{
	Kernel::M_TYPE mt = _basicMat;
	_maxNeigh  = third();
	_neighbors = new LatticeAtom *[_maxNeigh];
	for(unsigned i = 0; i < _maxNeigh; ++i)
		_neighbors[i] = 0;

	Mesh::LSNeiList neighbors;
	float minDist = _pDomain->_pLaPar[mt]->_neighborDistance[0];
	float second2 = _pDomain->_pLaPar[mt]->_neighborDistance[1] * _pDomain->_pLaPar[mt]->_neighborDistance[1];
	float first2 = minDist * minDist;
	minDist *= MIN_DIST_FACTOR;
	_pDomain->_pMesh->fillLatticeNeighbors(_coord, minDist, _pDomain->_pLaPar[mt]->_neighborDistance[2], neighbors);
	float minDist2 = first2 * MIN_DIST_FACTOR;

	bool good = true;
	for(Mesh::LSNeiList::iterator it=neighbors.begin(); it!=neighbors.end() && good; ++it)
	{
		LatticeAtom *itpLA = static_cast<LatticeAtom *>(it->_pLS);
		if(itpLA->getLatticeType() == getLatticeType() && itpLA != this)
		{
			int this2it, it2this, neighborhood;
			if(it->_dist2 < first2)
			{
				neighborhood = 0;
			}
			else if(it->_dist2 < second2)
			{
				neighborhood = 1;
			}
			else
			{
				neighborhood = 2;
			}
			this2it = itpLA->tryNeig(neighborhood, this, it->_dist2, minDist2);
			it2this = this->tryNeig(neighborhood, itpLA, it->_dist2, minDist2);
			if(this2it >= 0 && it2this >= 0) {
				itpLA->insertNeig(this, this2it, it->_dist2);
				this->insertNeig(itpLA, it2this, it->_dist2);
			}
			else if(this2it == NEIGH_NO_PLACE || it2this == NEIGH_NO_PLACE) {
				if(this2it == NEIGH_NO_PLACE) {
					itpLA->error(this, "neig");
				}
				if(it2this == NEIGH_NO_PLACE) {
					error(itpLA, "host");
				}
				int limit;
				if(neighborhood == 0) {
					limit = first();
				}
				else if(neighborhood == 1) {
					limit = second() - first();
				}
				else {
					limit = third() - second();
				}
				ERRORMSG("Inserting too many neighbors for neighborhood " << neighborhood + 1 << ", limit: " << limit);
			}
			else {} // Too close to an existing neighbour, skip silently.
		}
	}
}

void LatticeAtom::removeNeighbor(LatticeAtom *pLA)
{
	for(unsigned i=0; i<this->third(); ++i)
		if(_neighbors[i] == pLA)
		{
			_neighbors[i] = 0;
			return;
		}
	WARNINGMSG("Removing neighbor that does not exist!");
}

int LatticeAtom::tryNeig(unsigned in, LatticeAtom *pLA, float dist2, float minDist2)
{
	unsigned start = _maxNeigh;
	unsigned end = _maxNeigh;
	switch (in)
	{
	case 0:
		start = 0;
		end = this->first();
		break;
	case 1:
		start=this->first();
		end = this->second();
		break;
	case 2:
		start=this->second();
		end = this->third();
	}
	assert(start < end && end <= _maxNeigh);
	auto const pMesh = _pElement->getDomain()->_pMesh;
	for(unsigned i=0; i < _maxNeigh; ++i) {
		if(_neighbors[i] && pMesh->getDist2periodic(_neighbors[i]->_coord, pLA->_coord) < minDist2) {
			return NEIGH_TOO_CLOSE;
		}
	}
	for(unsigned i=start; i<end; ++i) {
		if(!_neighbors[i])
		{
			return i;
		}
	}
	return NEIGH_NO_PLACE;
}

void LatticeAtom::insertNeig(LatticeAtom * const pLA, int const where, float dist2) {
	_neighbors[where] = pLA;
	this->setNeiDist(where, dist2);
}

int LatticeAtom::getCoordination(unsigned i) const
{
	unsigned coordination = 0;
	unsigned start, end;
	switch (i)
	{
	case 0:
		start=0;
		end = this->first();
		break;
	case 1:
		start=this->first();
		end  =this->second();
		break;
	case 2:
		start=this->second();
		end  =this->third();
		break;
	default:
		start = end = 0;
		ERRORMSG("getCoordination() " << i);
		break;
	}
	for(unsigned i=start; i<end; ++i)
		if(_neighbors[i] && _neighbors[i]->getPerformed())
			coordination++;
	return coordination;
}

void LatticeAtom::displaceCoordinates(Kernel::Coordinates &c)
{
	static Kernel::Coordinates dummy(0,0,0);
	if(c == dummy)
		return;

	Kernel::MeshElement *pEle = _pElement;
	Kernel::Coordinates coord = _coord;
	if (_pDomain->_pMesh->jumpPosition(coord, c, pEle, 0, dummy, 0) != Kernel::Mesh::JUMP_OK)
		return;
	_coord = coord;
	if(pEle != _pElement)
	{
		_pDomain->_pMesh->remove(this);
		_pDomain->_pMesh->insert(this,pEle);
	}
	//update neighbors
	for(unsigned i=0; i < third(); ++i)
		if(_neighbors[i])
		{
			_neighbors[i]->removeNeighbor(this);
			_neighbors[i] = 0;
		}

	insertNeighbors();
}

void LatticeAtom::updateME4Epitaxy(Kernel::SubDomain *pSub)
{
	if(_pElement->getCrystallineLA() != 1)
		return; //not the first event performed in the cell.
	//first event. Check all the cells around and add more "stuff" if Gas is Found

	vector<Kernel::MeshElement *> neighbors;
	Kernel::M_TYPE mtN = Kernel::MAX_MATERIALS;
	_pDomain->_pMesh->fillAdjacentNeighbors(_pElement, neighbors, _pElement->getMaterial(), Kernel::MAX_MATERIALS);
	for(vector<Kernel::MeshElement *>::iterator itN = neighbors.begin(); itN != neighbors.end(); ++itN)
		if((*itN)->getFirstLS())
		{
			mtN = (*itN)->getFirstLS()->getBasicMat();
			break;
		}
	assert(mtN != Kernel::MAX_MATERIALS);

	for(vector<Kernel::MeshElement *>::iterator itN = neighbors.begin(); itN != neighbors.end(); ++itN)
		if( !(*itN)->getFirstLS() )
			_pDomain->_pLKMC->putLKMCAtoms(pSub, *itN, mtN, MOD_Epitaxy);
}

}
