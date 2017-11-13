/*
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

#include "MeshElement.h"
#include "Domain.h"
#include "SubDomain.h"
#include "Coordinates.h"
#include "Mesh.h"
#include "lkmc/Lattice.h"
#include "lkmc/LatticeAtom.h"
#include "lkmc/LKMCModel.h"
#include "RateManager.h"
#include "MeshParam.h"
#include "okmc/Interface.h"
#include "okmc/AlloyParam.h"
#include "okmc/MobileParticleParam.h"
#include "io/ParameterManager.h"

#include <cstring>

using namespace Kernel;
using std::vector;
using std::make_pair;
using std::pair;
using LKMC::LatticeAtom;

MeshElement::MeshElement(Domain *p) : _pDomain(p)
{
	_subDomainIdx = 0; //to be assigned later in Mesh
	_subDomainLevel = 0;
	_firstLS = 0;
	_firstPart = 0;
	_howManyParts = 0;
	_nonCrystallineLA = _crystallineLA = 0; //no LA atoms yet
	_AAtoms = 0;
    _BAtoms = 0;
	_stress = ublas::zero_vector<double>(6);
	_strain = ublas::zero_vector<double>(6);

	_V      = 0;
	_Eg		= 0;
	_index  = 0;
	_mat	= Kernel::MAX_MATERIALS;
	
	_amorphParts = 0;
}


MeshElement::~MeshElement()
{
}

void MeshElement::updateLAStatus(Kernel::SubDomain *pSub, LKMC::LKMCMode mode)
{
	_crystallineLA++;
	_nonCrystallineLA--;
	M_TYPE alloy = Domains::global()->PM()->getMaterial(getFirstLS()->getBasicMat())._alloy;
	if(_nonCrystallineLA == 0)
	{//element is completely crystalline
		M_TYPE toMat = getToMat(mode);
		_pDomain->_pMesh->changeMaterial(pSub, getIndex(), toMat);
		if(Domains::global()->PM()->isGasLike(_mat)) //epitaxy
		{
			LKMC::LatticeSite *pLS = getFirstLS();
			while(pLS)
			{
				if(pLS->getPType() == alloy)
                {
					decAAtoms();
					incBAtoms();
                }
				pLS = pLS->getNext();
			}
		}
	}
}

void MeshElement::updateUnLAStatus(Kernel::SubDomain *pSub, LKMC::LKMCMode mode)
{
	_crystallineLA--;
	_nonCrystallineLA++;
	if(_nonCrystallineLA == 1)
	{//first non-crystalline element in the element.
		M_TYPE toMat = getToMat(mode);
		assert(Domains::global()->PM()->isGasLike(toMat));
		_pDomain->_pMesh->changeMaterial(pSub, getIndex(), toMat);
		_BAtoms = 0;
	}
}

//when a new atoms shows up in a box, to check whether it is the first one, and then the toMat has to be updated...
M_TYPE MeshElement::getToMat(LKMC::LKMCMode mode)
{
	const IO::MaterialProperties &mat = Domains::global()->PM()->getMaterial(_mat);
	M_TYPE toMat = MAX_MATERIALS;

	switch(mode)
	{
	case LKMC::MOD_SPER:
		toMat = (mat._bAmorphous? mat._basicMaterial : mat._amorphMaterial);
		break;
	case LKMC::MOD_Epitaxy:
		if(Domains::global()->PM()->isGasLike(_mat))//it sounds like Epitaxy
		{
			vector<MeshElement *> neighbors;
			_pDomain->_pMesh->fillAdjacentNeighbors(this, neighbors, MAX_MATERIALS, MAX_MATERIALS);
			vector<MeshElement *>::iterator itN = neighbors.begin();
			for(; itN != neighbors.end(); ++itN)
			{
				if((*itN)->getFirstLS())
				{
					toMat = (*itN)->getFirstLS()->getBasicMat();
					break;
				}
			}
			if(itN == neighbors.end())
				ERRORMSG("Trying to update material, but not found for " << Domains::global()->PM()->getMaterialName(_mat));
		}
		else //etching, return gas
		{
			for(Kernel::M_TYPE mt = 0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
				if(Domains::global()->PM()->isGasLike(mt))
				{
					toMat = mt;
					break;
				}
		}

	}

	assert(toMat != MAX_MATERIALS);
	return toMat;
}

void MeshElement::insert(OKMC::Interface *pIF)
{
	_pDomain->_pMesh->trackLHF(_index, +1);
	_interfaces.push_back(pIF);
}

void MeshElement::remove(OKMC::Interface *pIF)
{
	_pDomain->_pMesh->trackLHF(_index, -1);
	for(std::vector<OKMC::Interface *>::iterator it = _interfaces.begin(); it!= _interfaces.end(); ++it)
		if(*it == pIF)
		{
			*it = _interfaces.back();
			_interfaces.pop_back();
			return;
		}
	ERRORMSG("Removing interface that was NOT inserted! Aing, aing!");
}

unsigned MeshElement::getHops(P_TYPE pt, unsigned st) const
{
	std::map<std::pair<P_TYPE, unsigned>, unsigned >::const_iterator it =
			_hops.find(make_pair(pt, st));
	if(it == _hops.end())
		return 0;
	return it->second;
}

Kernel::Coordinates MeshElement::getCoords(float rx, float ry, float rz) const
{
	Kernel::Coordinates m, M;
	_pDomain->_pMesh->getCorners(getIndex(), m, M);
	return Kernel::Coordinates(m._x + (M._x - m._x)*rx, m._y + (M._y - m._y)*ry, m._z + (M._z - m._z)*rz);
}

void MeshElement::getCorners(Kernel::Coordinates &m, Kernel::Coordinates &M) const
{
	_pDomain->_pMesh->getCorners(getIndex(), m, M);
}

double MeshElement::getVolume() const
{
	return _pDomain->_pMesh->getVolume(_index);
}

double MeshElement::getAlloyFraction() const {
	return _AAtoms ? double(_BAtoms) / double(getAtoms()) : 1.;
}

double MeshElement::getBasisFraction() const {
	return _BAtoms ? double(_AAtoms) / double(getAtoms()) : 1.;
}

double MeshElement::getEffectiveAlloyFraction() const {
    double wSides = _pDomain->_pAlPar->_smoothing[_mat][0];
    double wEdges = _pDomain->_pAlPar->_smoothing[_mat][1];
    double wCorners = _pDomain->_pAlPar->_smoothing[_mat][2];
    double wCentral = 1 - 6 * wSides - 12 * wEdges - 8 * wCorners;
    // Mirror Boundary conditions when different material or gas
    wCentral += (6 - _sides.size()) * wSides + (12 - _edges.size()) * wEdges + (8 - _corners.size()) * wCorners;

    double ret = wCentral * getAlloyFraction();

    for(vector<const MeshElement *>::const_iterator it = _sides.begin(); it != _sides.end(); ++it)
    	ret += wSides * (*it)->getAlloyFraction();
    for(vector<const MeshElement *>::const_iterator it = _edges.begin(); it != _edges.end(); ++it)
    	ret += wEdges * (*it)->getAlloyFraction();
    for(vector<const MeshElement *>::const_iterator it = _corners.begin(); it != _corners.end(); ++it)
    	ret += wCorners * (*it)->getAlloyFraction();
    return ret;
}

void MeshElement::decAAtoms() {
	if(!_AAtoms)
		ERRORMSG("Trying to remove an A atom from a cell but A atom counter is actually 0");
	_AAtoms--;
	if(!getAtoms())
	{
		Coordinates center;
		_pDomain->_pMesh->getCenter(center, _index);
		ERRORMSG("Both counters for the atoms in the cell (" << center << ") reached 0");
	}
}

void MeshElement::decBAtoms() {
	if(!_BAtoms)
		ERRORMSG("Trying to remove an B atom from a cell but B atom counter is actually 0");
	_BAtoms--;
	if(!getAtoms())
	{
		Coordinates center;
		_pDomain->_pMesh->getCenter(center, _index);
		ERRORMSG("Both counters for the atoms in the cell (" << center << ") reached 0");
	}
}

void MeshElement::amorphizeME(Kernel::SubDomain *pSub)
{
	M_TYPE toMat = getToMat(LKMC::MOD_SPER);
	_pDomain->_pMesh->changeMaterial(pSub,getIndex(), toMat);
}

void MeshElement::restart(std::istream &is)
{
	is >> _BAtoms >> _AAtoms;
}

void MeshElement::restart(std::ostream &os) const
{
	os << _BAtoms << " " << _AAtoms << " ";
}

