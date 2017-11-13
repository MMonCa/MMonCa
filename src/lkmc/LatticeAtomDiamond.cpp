/*
 * LatticeAtomDiamond.cpp
 *
 *  Created on: Oct 10, 2011
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

#include "LatticeAtomDiamond.h"

#include "kernel/Mesh.h"
#include "kernel/RateManager.h"
#include "kernel/SubDomain.h"
#include "domains/Global.h"
#include "io/ParameterManager.h"
#include "LatticeDiamondParam.h"
#include "LatticeAtomDSPER.h"
#include "LatticeAtomDEpi.h"
#include "LatticeAtomD2Epi.h"

using Kernel::Mesh;
using Kernel::Coordinates;
using std::vector;
using std::pair;
using std::map;

namespace LKMC {

LatticeAtomDiamond::LatticeAtomDiamond(LASTATE per, Kernel::Domain *p, Kernel::M_TYPE basicMt, Kernel::P_TYPE type, const Coordinates &c, unsigned orientation) :
		LatticeAtom(per, p, basicMt, type, c), _orientation(orientation)
{
}

LatticeAtomDiamond::LatticeAtomDiamond(std::istream &is) : LatticeAtom(is)
{
	is >> _orientation;
}

void LatticeAtomDiamond::restart(std::ostream &os) const
{
	LatticeAtom::restart(os);
	os << _orientation << " ";
}

//Conditions for "isAllowed" are:
// there is not another particle around in a reasonable radius
// the particle has not many first neighbors (more than 5?)
bool LatticeAtomDiamond::isAllowed(const Coordinates &final, const Kernel::MeshElement *pEle, LKMCMode mode) const
{
	const float neiDist0 = _pDomain->_pLaPar[_basicMat]->_neighborDistance[0];
	const static float maxDist2 = 0.20*0.20;
	IO::ParameterManager *pPM = Domains::global()->PM();
	Kernel::M_TYPE mt = pEle->getMaterial();	
	if( pPM->getMaterial(mt)._basicMaterial != _basicMat &&
	    ( (mode == MOD_Epitaxy && !pPM->isGasLike(mt)) || mode == MOD_SPER))
	    	return false;
	Mesh::LSNeiList neiList;
	_pDomain->_pMesh->fillLatticeNeighbors(final, neiDist0, neiList);
	if(neiList.size() > 5)
		return false;
	for(Mesh::LSNeiList::iterator it=neiList.begin(); it != neiList.end(); ++it)
		if(it->_dist2 < maxDist2)
			return false;
	return true;
}

bool LatticeAtomDiamond::isDefective() const
{
	unsigned howMany = 0;
	for(unsigned i=0; i<first(); ++i)
		if(_neighbors[i])
			howMany++;
	return howMany > 4;
}

//helper to create atoms from a neighboring list
template <typename AtomType>
void LatticeAtomDiamond::createFromNeiList(Kernel::SubDomain *pSub, LASTATE st, const vector<pair<Coordinates, unsigned> > &neis, AtomType *pAmo, LKMCMode mode, map<int, LatticeAtomDiamond *> &toUpdate) const
{
	Coordinates dummy(0,0,0);
	for(vector<pair<Coordinates, unsigned> >::const_iterator it=neis.begin(); it != neis.end(); ++it)
	{
		Mesh::LSNeiList neiList;
		Kernel::MeshElement *pEle = pAmo->_pElement;
		Kernel::Coordinates final = pAmo->_coord;
		if(_pDomain->_pMesh->jumpPosition(final, it->first, pEle, 0, dummy, 0) == Mesh::JUMP_OK && isAllowed(final, pEle, mode))
		{
			pEle->incLA(st == LS_PERFORMED);
			Kernel::M_TYPE alloy = Domains::global()->PM()->getMaterial(_basicMat)._alloy;
			Kernel::P_TYPE myPt = _type;
			if(mode == MOD_SPER && alloy != Kernel::MAX_MATERIALS) //alloy is defined
			{
				float prob = pEle->getBAtoms() / float(pEle->getAtoms());
				myPt =  (pSub->_rng.rand() < prob? alloy : _type);
			}
			AtomType *pLA = new AtomType(st, _pDomain, _basicMat, myPt, final, it->second);  //orientation is a "guess" for am. atoms
			_pDomain->_pRM->insert(pLA, pEle);
			toBeUpdated(pLA, toUpdate);
		}
	}
}

//the list of atoms to be updated are:
// the atom
// the first neighbors of the atom
// the second neighbors of every first neighbor. (note this is bigger than the third neighbors of the atom)
// the first neighbors of the second neighbors
// _number is used here to avoid introducing a dependency on the order of the pointers
void LatticeAtomDiamond::toBeUpdated(LatticeAtomDiamond *p, map<int, LatticeAtomDiamond *> &toUpdate)
{
	toUpdate[p->_number] = p;
	for(unsigned i=0; i<p->second(); ++i)
	{
		LatticeAtomDiamond *nei = static_cast<LatticeAtomDiamond *>(p->_neighbors[i]);
		if(nei)
		{
			toUpdate[nei->_number] = nei;
			unsigned end = (i< p->first() ? p->second() : p->first());
			for(unsigned j=0; j<end; ++j)
			{
				LatticeAtomDiamond *neinei = static_cast<LatticeAtomDiamond *>(nei->_neighbors[j]);
				if(neinei)
					toUpdate[neinei->_number] = neinei;
			}
		}
	}
}

//explicit instantiations
template void LatticeAtomDiamond::createFromNeiList(Kernel::SubDomain *, LASTATE, const vector<pair<Coordinates, unsigned> > &, LatticeAtomDSPER *, LKMCMode, map<int, LatticeAtomDiamond *> &) const;
template void LatticeAtomDiamond::createFromNeiList(Kernel::SubDomain *, LASTATE, const vector<pair<Coordinates, unsigned> > &, LatticeAtomDEpi  *, LKMCMode, map<int, LatticeAtomDiamond *> &) const;
template void LatticeAtomDiamond::createFromNeiList(Kernel::SubDomain *, LASTATE, const vector<pair<Coordinates, unsigned> > &, LatticeAtomD2Epi *, LKMCMode, map<int, LatticeAtomDiamond *> &) const;
}
