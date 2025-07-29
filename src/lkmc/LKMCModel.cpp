/*
 * LKMCModel.cpp
 *
 *  Created on: May 11, 2011
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

#include "LKMCModel.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "kernel/RateManager.h"
#include "kernel/MeshParam.h"
#include "kernel/SubDomain.h"
#include "io/Diagnostic.h"
#include "LatticeAtomDSPER.h"
#include "LatticeAtomDEpi.h"
#include "LatticeAtomD2Epi.h"
#include "LatticeAtomBCC1.h"
#include "LatticeParam.h"
#include "Lattice.h"
#include "io/ParameterManager.h"

using Kernel::Coordinates;
using std::vector;
using std::pair;
using Kernel::MeshElement;
using Kernel::Mesh;
using std::map;

namespace LKMC {

LKMCModel::LKMCModel(bool bFromStream, Kernel::Domain *pDomain)
{
	_pDomain = pDomain;
	_pDomain->_pLKMC = this;
	if(!bFromStream) //otherwise the restart process re-creates the atoms.
	{
		epitaxy();
		sper();
	}
}

//called to fill an empty box.
unsigned LKMCModel::putLKMCAtoms(Kernel::SubDomain *pSub, MeshElement *pEle, Kernel::M_TYPE basicMat, LKMCMode mode) const
{
	unsigned na = 0;
	IO::MaterialProperties mprop = Domains::global()->PM()->getMaterial(basicMat);
	Coordinates m, M;
	_pDomain->_pMesh->getCorners(pEle->getIndex(), m, M);
	vector<Lattice::LatticeInformation> atoms;
	_pDomain->_pLat[basicMat]->fill(m, M, atoms, false);

	LatticeAtom::LASTATE state = LatticeAtom::LS_AVAILABLE;
	if(mode == MOD_SPER)
		state = (Domains::global()->PM()->getMaterial(pEle->getMaterial())._bAmorphous == false ? LatticeAtom::LS_PERFORMED : LatticeAtom::LS_AVAILABLE);
	else
		state = (Domains::global()->PM()->getMaterial(pEle->getMaterial())._bEpitaxy? LatticeAtom::LS_PERFORMED : LatticeAtom::LS_AVAILABLE);

	for(vector<Lattice::LatticeInformation>::iterator it = atoms.begin(); it!=atoms.end();++it)
	{
		LatticeSite *pLS = _pDomain->_pMesh->findLS(it->_coords);
		Kernel::P_TYPE pt = (mprop._binary ? mprop._pt[it->_type] : mprop._pt[0]);
		if(pLS == 0)
		{
			LatticeAtom * pLA = 0;
			na++;
			if(mode == MOD_SPER)
			{
				if(_pDomain->_pLaPar[basicMat]->_type == LatticeParam::DIAMOND)
					pLA = new LKMC::LatticeAtomDSPER(state, _pDomain, basicMat, pt, it->_coords, it->_orientation);
				else
					ERRORMSG("Atom type not allowed for SPER");
			}
			else
			{
				if(_pDomain->_pLaPar[basicMat]->_type == LatticeParam::DIAMOND)
					pLA = new LKMC::LatticeAtomDEpi(state, _pDomain, basicMat, pt, it->_coords, it->_orientation);
				else if(_pDomain->_pLaPar[basicMat]->_type == LatticeParam::DIAMOND2)
					pLA = new LKMC::LatticeAtomD2Epi(state, _pDomain, basicMat, pt, it->_coords, it->_orientation);
				else if(_pDomain->_pLaPar[basicMat]->_type == LatticeParam::BCC)
					pLA = new LKMC::LatticeAtomBCC1(state, _pDomain, basicMat, pt, it->_coords);
				else
					ERRORMSG("LKMCModel:Lattice Type NOT DEFINED for LKMC");
			}
			pLA->getElement()->incLA(state == LatticeAtom::LS_PERFORMED);
			_pDomain->_pRM->insert(pLA, pLA->getElement());
		}
		else //atom does not need to be in exactly the same element...
		{
			LatticeAtom *pLA = static_cast<LatticeAtom *>(pLS);
			pLA->setState(state);
			pEle->getDomain()->_pRM->update(pLA, pLA->getElement());
		}
	}
	return na;
}

//Cleans lattice atoms as the interface changes
void LKMCModel::cleanLKMCAtoms(MeshElement *pEle, LKMCMode mode) const
{
	const IO::ParameterManager *pPM = Domains::global()->PM();
	vector<MeshElement *> toRemove;

	if(mode == MOD_SPER)
	{
		Kernel::M_TYPE basicMat = pPM->getMaterial(pEle->getMaterial())._basicMaterial;
		vector<MeshElement *> elements;
		pEle->getDomain()->_pMesh->fillAdjacentNeighbors(pEle, elements, Kernel::MAX_MATERIALS, basicMat);
		elements.push_back(pEle);
		//For each neighbor having ALL of its neighbors with the "same crystallinity", remove the lattice
		for(vector<MeshElement *>::iterator it=elements.begin(); it!=elements.end(); ++it)
		{
			//get the neighboring neighbors
			vector<MeshElement *> neiEles;
			(*it)->getDomain()->_pMesh->fillAdjacentNeighbors(*it, neiEles, Kernel::MAX_MATERIALS, basicMat);
			bool bRemove = true;
			//if just one has "different crystallinity", do not remove anything.
			bool bAmor = pPM->isAmorphous((*it)->getMaterial());
			for(vector<MeshElement *>::iterator neiIt=neiEles.begin(); neiIt!=neiEles.end(); ++neiIt)
				if( (bAmor == false && (*neiIt)->getMaterial() != (*it)->getMaterial()) ||
					(bAmor == true && ((*it)->getCrystallineLA() != 0 || (*neiIt)->getCrystallineLA() != 0 )))
				{
					bRemove = false;
					break;
				}
			if(bRemove)
				toRemove.push_back(*it);
		}
	}
	if(mode == MOD_Epitaxy)
	{
		vector<MeshElement *> elements;
		pEle->getDomain()->_pMesh->fillAdjacentNeighbors(pEle, elements, Kernel::MAX_MATERIALS, Kernel::MAX_MATERIALS);
		Kernel::M_TYPE mt = pEle->getMaterial();
		if(pPM->isGasLike(mt))
			mt = pEle->getToMat(MOD_Epitaxy);
		assert(!pPM->isGasLike(mt));
		elements.push_back(pEle);
		//For each neighbor having ALL of its neighbors the same, remove the lattice
		for(vector<MeshElement *>::iterator it=elements.begin(); it!=elements.end(); ++it)
		{
			if((*it)->getCrystallineLA() && (*it)->getNonCrystallineLA())
				continue;
			//get the neighboring neighbors
			vector<MeshElement *> neiEles;
			(*it)->getDomain()->_pMesh->fillAdjacentNeighbors(*it, neiEles, Kernel::MAX_MATERIALS, Kernel::MAX_MATERIALS);
			bool bRemove = true;
			//if just one has "different crystallinity", do not remove anything.
			for(vector<MeshElement *>::iterator neiIt=neiEles.begin(); neiIt!=neiEles.end(); ++neiIt)
			{
				if((*neiIt)->getCrystallineLA() && (*neiIt)->getNonCrystallineLA()) // half recrystallized
				{
					bRemove = false;
					break;
				}
				if(((*it)->getCrystallineLA() == 0 && (*neiIt)->getCrystallineLA() != 0) ||
				   ((*it)->getNonCrystallineLA() == 0 && (*neiIt)->getNonCrystallineLA() != 0))
				{
					bRemove = false;
					break;
				}
			}
			if(bRemove)
				toRemove.push_back(*it);
		}
	}

	for(vector<MeshElement *>::iterator it=toRemove.begin(); it!=toRemove.end(); ++it) //proceed to remove stuff...
	{
		vector<LatticeSite *> atoms;
		LatticeSite *pLS = (*it)->getFirstLS();
		while(pLS)
		{
			atoms.push_back(pLS);
			pLS = pLS->getNext();
		}
		for(vector<LatticeSite *>::iterator laIt = atoms.begin(); laIt!=atoms.end(); ++laIt)
		{
			LatticeAtom *pLA = static_cast<LatticeAtom *>(*laIt);
			pLA->getElement()->decLA(pLA->getPerformed());
			_pDomain->_pRM->remove(pLA, pLA->getElement());
			delete pLA;
		}
	}
}

void LKMCModel::epitaxy() const
{
	Mesh * pMesh = _pDomain->_pMesh;
	const IO::ParameterManager *pPM = Domains::global()->PM();

	for(Mesh::iterator it=pMesh->begin(); it!=pMesh->end(); ++it)
	{
		Kernel::SubDomain *pSub = _pDomain->_pRM->getSubDomain(it->getSubDomainIdx());
		Kernel::M_TYPE mt = it->getMaterial();
		if(pPM->isGasLike(mt))  //Epitaxy
		{
			vector<MeshElement *> neighbors;
			pMesh->fillAdjacentNeighbors(&*it, neighbors, Kernel::MAX_MATERIALS, Kernel::MAX_MATERIALS);
			for(vector<MeshElement *>::iterator itN = neighbors.begin(); itN != neighbors.end(); ++itN)
			{
				Kernel::M_TYPE mtN = (*itN)->getMaterial();
				if(pPM->getMaterial(mtN)._bEpitaxy)
				{
					putLKMCAtoms(pSub, *itN, mtN, MOD_Epitaxy);
					putLKMCAtoms(pSub, &*it, mtN, MOD_Epitaxy);
				}
			}
		}
	}
	//Refinement
	for(Mesh::iterator it=pMesh->begin(); it!=pMesh->end(); ++it)
	{
		LatticeSite *pLS = it->getFirstLS();
		while(pLS)
		{
			LatticeAtomDEpi * pLA = dynamic_cast<LatticeAtomDEpi *>(pLS);
			if (pLA && !pLA->getPerformed() && pLA->getCoordination(0) == 0) //not at the interface
			{
				pLA->getElement()->decLA(pLA->getPerformed());
				_pDomain->_pRM->remove(pLA, pLA->getElement());
				pLS = pLS->getNext();
				delete pLA;
			}
			else {
				pLS = pLS->getNext();
			}
		}
	}
}

void LKMCModel::sper() const
{
	LOWMSG2("Checking SPER...");
	std::vector <Kernel::MeshElement *> elems;
	for(Mesh::iterator it=_pDomain->_pMesh->begin(); it!=_pDomain->_pMesh->end(); ++it)
		if(it->getInterfaces().size())
			elems.push_back(&*it);
	localSPER(_pDomain->_pRM->getSubDomain(0), elems);
	LOWMSG(" Done.");
}

//this function is only called when crystalline cells are amorphized.
//if the material is already amorphous, and it is called to amorphize it again, its behavior is unpredicted
void LKMCModel::localSPER(Kernel::SubDomain *pSub, const std::vector<Kernel::MeshElement *> &elems) const
{
	std::map<unsigned, Kernel::MeshElement *> cells,cellsToUpdate;
	const IO::ParameterManager *pPM = Domains::global()->PM();
	Mesh * pMesh = _pDomain->_pMesh;
	for(std::vector<Kernel::MeshElement *>::const_iterator it=elems.begin(); it!=elems.end(); ++it)
	{
		Kernel::M_TYPE mt = (*it)->getMaterial();
		if(pPM->getMaterial(mt)._bAmorphous)  //sper
		{
			cellsToUpdate[(*it)->getIndex()] = *it;
			cells[(*it)->getIndex()] = *it;
			Kernel::M_TYPE basMat = pPM->getMaterial(mt)._basicMaterial;
			vector<MeshElement *> neighbors;
			pMesh->fillAdjacentNeighbors(*it, neighbors, Kernel::MAX_MATERIALS, Kernel::MAX_MATERIALS);
			for(vector<MeshElement *>::iterator itN = neighbors.begin(); itN != neighbors.end(); ++itN)
			{
				Kernel::M_TYPE mtN = (*itN)->getMaterial();
				if(pPM->getMaterial(mtN)._basicMaterial == basMat)
					cellsToUpdate[(*itN)->getIndex()] = *itN;
				if(!pPM->getMaterial(mtN)._bAmorphous && pPM->getMaterial(mtN)._basicMaterial == basMat)
					if((*itN)->getFirstLS() == 0 || static_cast<LatticeAtom *>((*itN)->getFirstLS())->getPerformed() == false)
						cells[(*itN)->getIndex()] = *itN;
			}
		}
	}
	//creation
	unsigned na = 0;
	for(std::map<unsigned, Kernel::MeshElement *>::iterator it=cells.begin(); it!=cells.end(); ++it)
		na += putLKMCAtoms(pSub, it->second, pPM->getMaterial(it->second->getMaterial())._basicMaterial, MOD_SPER);
	LOWMSG2("  +" << na);
	//Refinement, it removes atoms!
	na = 0u;
	for(std::map<unsigned, Kernel::MeshElement *>::const_iterator it=cellsToUpdate.begin(); it!=cellsToUpdate.end(); ++it)
	{
		LatticeSite *pLS = it->second->getFirstLS();
		while(pLS)
		{
			LatticeAtomDSPER * pLA = dynamic_cast<LatticeAtomDSPER *>(pLS);
			if (pLA && !pLA->getPerformed() && pLA->getCoordination(0) == 0) //not at the interface
			{
				pLA->getElement()->decLA(pLA->getPerformed());
				_pDomain->_pRM->remove(pLA, pLA->getElement());
				pLS = pLS->getNext();
				delete pLA;
				na++;
			}
			else {
				pLS = pLS->getNext();
			}
		}
	}
	LOWMSG2(" / -" << na << " atoms.");
	//build updated list now that atoms were removed.
	std::map<int, LatticeAtomDiamond *> toUpdate;
	for(std::map<unsigned, Kernel::MeshElement *>::iterator it=cellsToUpdate.begin(); it!=cellsToUpdate.end(); ++it)
	{
		LatticeSite *pLS = it->second->getFirstLS();
		while(pLS)
		{
			LatticeAtomDSPER * pLA = dynamic_cast<LatticeAtomDSPER *>(pLS);
			if(pLA)
				LatticeAtomDiamond::toBeUpdated(pLA, toUpdate);
			pLS = pLS->getNext();
		}
	}
	//update
	for(std::map<int, LatticeAtomDiamond *>::iterator it=toUpdate.begin(); it!=toUpdate.end(); ++it)
		it->second->getDomain()->_pRM->update(it->second, it->second->getElement());
}

}
