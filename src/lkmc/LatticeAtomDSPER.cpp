/*
 * LatticeAtomDSPER.cpp
 *
 *  Created on: Jul 15, 2014
 *      Author: ignacio.martin@imdea.org
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

#include "LatticeAtomDSPER.h"
#include "LatticeDiamondParam.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/RateManager.h"
#include "kernel/MeshElement.h"
#include "kernel/Mesh.h"
#include "LatticeDiamond.h"
#include "io/ParameterManager.h"
#include "kernel/Constants.h"

using std::map;
using std::vector;
using Kernel::Coordinates;
using std::pair;
using Kernel::Mesh;

namespace LKMC {

LatticeAtomDSPER::LatticeAtomDSPER(LASTATE st, Kernel::Domain *p, Kernel::M_TYPE basicMt, Kernel::P_TYPE type, const Kernel::Coordinates &c, unsigned orientation)
 : LatticeAtomDiamond(st, p, basicMt, type, c, orientation)
{
	_plane = LatticeDiamondParam::PNONE;
	insertNeighbors();  //this function cannot be in the base class because it needs the virtual classes working
}

LatticeAtomDSPER::LatticeAtomDSPER(std::istream &is) : LatticeAtomDiamond(is)
{
	_plane = LatticeDiamondParam::PNONE;
	insertNeighbors();
}

void LatticeAtomDSPER::restart(std::ostream &os) const
{
	LatticeAtomDiamond::restart(os);
}

void LatticeAtomDSPER::perform(Kernel::SubDomain *pSub, unsigned eventType)
{
	if(eventType != 0)
		ERRORMSG("Unknown LKMC event");
	pSub->_evLog.performed(0, Event::LATTICEATOM, LSDT_DIAMOND_SPER, 0, _plane, 0);
	_pDomain->_pRM->setDepthLA(_coord._x);

	map<int, LatticeAtomDiamond *> toUpdate;
	switch(_plane)
	{
	case LatticeDiamondParam::P111:
		event111(pSub, toUpdate);
		break;
	case  LatticeDiamondParam::P110:
		event110(pSub, toUpdate);
		break;
	case LatticeDiamondParam::P100_10:
	case LatticeDiamondParam::P100_9:
	case LatticeDiamondParam::P100_8:
	case LatticeDiamondParam::P100_7:
	case LatticeDiamondParam::P100_6: //100
		event100(pSub, toUpdate);
		break;
	default:
		ERRORMSG("Strange event");
		break;
	}

	for(map<int, LatticeAtomDiamond *>::iterator it=toUpdate.begin(); it!=toUpdate.end(); ++it)
		_pDomain->_pRM->update(it->second, static_cast<LatticeAtomDSPER *>(it->second)->_pElement);
	if(_pElement->getNonCrystallineLA() == 0)
		_pDomain->_pLKMC->cleanLKMCAtoms(_pElement, MOD_SPER);
}

void LatticeAtomDSPER::event100(Kernel::SubDomain *pSub, map<int, LatticeAtomDiamond *> &toUpdate)
{
	Kernel::M_TYPE mt = _pElement->getMaterial();
	const LatticeDiamond *pLat = static_cast<const LatticeDiamond *>(_pDomain->_pLat[mt]);
	vector<pair<Coordinates, unsigned> > neis;
	//obtain director vectors and neighbors
	Coordinates v[2];
	int nC=0;
	for(unsigned i=0; i<FIRSTN; ++i)
		if(_neighbors[i] && _neighbors[i]->getPerformed())
		{
			v[nC] = _neighbors[i]->getCoordinates();
			_pDomain->_pMesh->setPeriodicRelative(_coord, v[nC]);
			nC++;
			if(nC == 2)
				break;
		}
	_orientation = pLat->getNeis2(pSub, 0, v[0], v[1], neis);

	createFromNeiList(pSub, LS_AVAILABLE, neis, this, MOD_SPER, toUpdate);
	_state = LS_PERFORMED;
	_pElement->updateLAStatus(pSub, MOD_SPER);
	toBeUpdated(this, toUpdate);
}

void LatticeAtomDSPER::event110(Kernel::SubDomain *pSub, map<int, LatticeAtomDiamond *> &toUpdate)
{
	Kernel::M_TYPE mt = _pElement->getMaterial();
	const LatticeDiamond *pLat = static_cast<const LatticeDiamond *>(_pDomain->_pLat[mt]);
	vector<LatticeAtomDSPER *> chains110, crys110;
	LatticeAtomDSPER *cry[2] = { 0, 0 };
	unsigned orient = 15;
	for(int i=0; i<FIRSTN; ++i)
	{
		LatticeAtomDSPER *pNei = static_cast<LatticeAtomDSPER *>(_neighbors[i]);
		if(pNei && pNei->getPerformed()) //if "crystalline" neighbor
		{
			cry[0] = pNei;
			orient = pNei->_orientation;
		}
	}
	assert(cry[0]);
	for(int i=0; i<FIRSTN; ++i)
	{
		LatticeAtomDSPER *pNei = static_cast<LatticeAtomDSPER *>(_neighbors[i]);
		if(pNei && !pNei->getPerformed())
			for(int j=0; j<FIRSTN; ++j)
			{
				LatticeAtomDSPER *pNeiNei = static_cast<LatticeAtomDSPER *>(pNei->_neighbors[j]);
				if(pNeiNei && pNeiNei->getPerformed() && pNeiNei->_orientation/2 == orient/2)
				{
					assert(pNeiNei != this);
					chains110.push_back(pNei);
					crys110.push_back(pNeiNei);
					break;
				}
			}
	}
	assert(chains110.size());
	const unsigned whichOne = unsigned(chains110.size()*pSub->_rng.rand());
	LatticeAtomDSPER * pDef = chains110[whichOne];
	cry[1] = crys110[whichOne];

	//This amorphous atom with its crystalline neighbor
	Coordinates v[2] = { cry[0]->getCoordinates(), pDef->getCoordinates() };
	_pDomain->_pMesh->setPeriodicRelative(_coord, v[0]);
	_pDomain->_pMesh->setPeriodicRelative(_coord, v[1]);
	vector<pair<Coordinates, unsigned> > neis;
	_orientation = pLat->getNeis2(pSub, 0, v[0], v[1], neis);
	createFromNeiList(pSub, LS_AVAILABLE, neis, this, MOD_SPER, toUpdate);

	//The second amorphous atom with its crystalline neighbor
	neis.clear();
	v[0] = _coord;
	v[1] = cry[1]->getCoordinates();
	_pDomain->_pMesh->setPeriodicRelative(pDef->getCoordinates(), v[0]);
	_pDomain->_pMesh->setPeriodicRelative(pDef->getCoordinates(), v[1]);
	pDef->_orientation = pLat->getNeis2(pSub, 0, v[0], v[1], neis);
	createFromNeiList(pSub, LS_AVAILABLE, neis, pDef, MOD_SPER, toUpdate);

	//not updating the Poisson nodes now...
	_state = pDef->_state = LS_PERFORMED;
	_pElement->updateLAStatus(pSub, MOD_SPER);
	pDef->_pElement->updateLAStatus(pSub, MOD_SPER);

	toBeUpdated(this, toUpdate);
	toBeUpdated(pDef, toUpdate);
}

void LatticeAtomDSPER::event111(Kernel::SubDomain *pSub, map<int, LatticeAtomDiamond *> &toUpdate)
{
	Kernel::M_TYPE mt = _pElement->getMaterial();
	const LatticeDiamond *pLat = static_cast<const LatticeDiamond *>(_pDomain->_pLat[mt]);
	//First obtain the two atoms that are crystalline
	LatticeAtomDSPER *cry[2] = {0, 0};
	for(int i=0; i<FIRSTN; ++i)
		if(_neighbors[i] && _neighbors[i]->getPerformed())
		{
			cry[0] = static_cast<LatticeAtomDSPER *>(_neighbors[i]);
			break;
		}

	vector<LatticeAtomDSPER *> candidatesCrys, candidatesAmor;
	for(int i=FIRSTN; i<SECONDN; ++i) //because the first neighbor in the 111 chain might not exist
	{
		LatticeAtomDSPER *pNei = static_cast<LatticeAtomDSPER *>(_neighbors[i]);
		if(pNei && !pNei->getPerformed()) //if "amorphous" neighbor
		{
			//unsigned coord=0;
			for(int j=0; j<FIRSTN; ++j)
				if(pNei->_neighbors[j] && pNei->_neighbors[j]->getPerformed())
				{
					cry[1] =  static_cast<LatticeAtomDSPER *>(pNei->_neighbors[j]);
					//coord++;
				}
			if(/*coord == 1 && */cry[1] && std::find(_neighbors+SECONDN, _neighbors+THIRDN, cry[1]) != _neighbors+THIRDN)
			{
				candidatesCrys.push_back(cry[1]);			//be sure it forms an "hexagon", i.e., it is a third neighbor.
				candidatesAmor.push_back(pNei);
			}
		}
	}
	if(candidatesCrys.empty())
		return;
	unsigned idx = unsigned(candidatesCrys.size()*pSub->_rng.rand());
	cry[1] = candidatesCrys[idx];
	//Create all atoms in the chain, including possible twins.
	LatticeAtomDSPER *amo[3] = {this, candidatesAmor[idx], 0};

	//check if there is already an amo[2]
	for(unsigned m=0; m<FIRSTN; ++m)
		for(unsigned j=0; j<FIRSTN; ++j)
			if(amo[0]->_neighbors[m] && amo[0]->_neighbors[m] == amo[1]->_neighbors[j])
			{
				amo[2] = static_cast<LatticeAtomDSPER *>(amo[0]->_neighbors[m]);
				if(amo[2]->getPerformed() == true)
					WARNINGMSG("Crystalline neighbor found for a <111> recrystallization!");
				//delete it. It might be in a wrong position.
				amo[2]->getElement()->decLA(amo[2]->getPerformed());
				_pDomain->_pRM->remove(amo[2], amo[2]->getElement());
				delete amo[2];
				break;
			}
	//create amo[2]
	{
		amo[2] = 0;
		Coordinates v[2] = { cry[0]->getCoordinates(), cry[1]->getCoordinates() };
		_pDomain->_pMesh->setPeriodicRelative(amo[0]->getCoordinates(), v[0]);
		_pDomain->_pMesh->setPeriodicRelative(amo[1]->getCoordinates(), v[1]);
		vector<pair<Coordinates, unsigned> > neis;
		int orientation = pLat->getNeis1(pSub, 0, v[0], v[1], neis, cry[0]->_orientation, cry[1]->_orientation); //might produce twin defects!
		createFromNeiList(pSub, LS_AVAILABLE, neis, amo[0], MOD_SPER, toUpdate);
		//check for amo[2]
		for(unsigned m=0; m<FIRSTN; ++m)
			for(unsigned j=0; j<FIRSTN; ++j)
				if(amo[0]->_neighbors[m] && amo[0]->_neighbors[m]== amo[1]->_neighbors[j])
				{
					amo[2] = static_cast<LatticeAtomDSPER *>(amo[0]->_neighbors[m]);
					amo[2]->_orientation = orientation;
					assert(!amo[2]->getPerformed());
				}

		//now for amo[0] and amo[1]
		for(unsigned i=0; i<2; ++i)
		{
			neis.clear();
			Coordinates v[2]; //temporal storage of coordinates needed to ask getNeis;
			v[0] = cry[i]->getCoordinates();
			_pDomain->_pMesh->setPeriodicRelative(amo[i]->getCoordinates(), v[0]);
			if(amo[2]) //obtain new neighboring positions
			{
				v[1] = amo[2]->getCoordinates();
				_pDomain->_pMesh->setPeriodicRelative(amo[i]->getCoordinates(), v[1]);
				amo[i]->_orientation = pLat->getNeis2(pSub, 0, v[0], v[1], neis);
				//for all neighboring positions try to create the atom if it is not there.
				createFromNeiList(pSub, LS_AVAILABLE, neis, amo[i], MOD_SPER, toUpdate);
			}
		}
	}
	if(amo[2]) //due to defects it might not have been created
	{
		vector<pair<Coordinates, unsigned> > neis;
		Coordinates v[2] = { amo[0]->getCoordinates(), amo[1]->getCoordinates() };
		_pDomain->_pMesh->setPeriodicRelative(amo[2]->getCoordinates(), v[0]);
		_pDomain->_pMesh->setPeriodicRelative(amo[2]->getCoordinates(), v[1]);
		pLat->getNeis2(pSub, 0, v[0], v[1], neis);
		createFromNeiList(pSub, LS_AVAILABLE, neis, amo[2], MOD_SPER, toUpdate);
	}

	//Start recrystallizing atoms and update everything.
	for(int j=0; j<3; ++j)
	{
		if(!amo[j])
			continue;
		amo[j]->_state = LS_PERFORMED;
		amo[j]->_pElement->updateLAStatus(pSub, MOD_SPER);
		toBeUpdated(amo[j], toUpdate);
	}
}

float LatticeAtomDSPER::getRate(unsigned eventType, float kT) const
{
	if (getPerformed() || eventType != 0)
		return 0;

	const LatticeDiamondParam * param = static_cast<const LatticeDiamondParam *>(_pDomain->_pLaPar[_basicMat]);
	unsigned plane = LatticeDiamondParam::PNONE;
	unsigned coord = 0, orient1=0, orient2=0;
	Kernel::ID neighborhood;
	neighborhood._mt = _basicMat;
	for(unsigned i=0; i<FIRSTN; ++i)
	{
		LatticeAtomDSPER *pLAD = static_cast<LatticeAtomDSPER *>(_neighbors[i]);
		if(pLAD && pLAD->getPerformed())
		{
			coord++;
			if(coord == 1)
				orient1 = orient2 = pLAD->_orientation/2;
			if(orient1 != pLAD->_orientation/2)
				orient2 = pLAD->_orientation/2;

			Kernel::P_TYPE pt = _neighbors[i]->getPType();
			std::map<Kernel::P_TYPE, unsigned>::iterator it = neighborhood._pt.find(pt);
			if(it == neighborhood._pt.end())
				neighborhood._pt[pt] = 1;
			else
				it->second++;
		}
	}
	plane = getPlane(coord);
	if(plane == LatticeDiamondParam::PNONE)
		return 0;
	float ener = 0;
	//compute energy.
	if(neighborhood._pt.size() != 0) //requires at least one atom to attach to.
	{
		unsigned hash = param->_sper[_type]._enerSPERFirst.cluster2hash(neighborhood);
		if(hash == param->_sper[_type]._enerSPERFirst.invalidHash())
			WARNINGMSG("Diamond SPER: Parameters not specified for " <<
				Domains::global()->PM()->getParticleName(neighborhood._mt, _type) << "+" <<
				Domains::global()->PM()->getEpiIDName(neighborhood) << ". Taking 0.");
		else
			ener += param->_sper[_type]._enerSPERFirst._map[hash];
	}
	ener += _pElement->strain_xy()*_pElement->strain_xy()*param->_shearEffect;
	if(orient1 == orient2)
		ener -= stressCorrection_2ndOrder(plane, coord, kT);
	_plane = plane;
	double electrostaticFactor = 1.;
	if (param->_electrostaticModel == LatticeDiamondParam::GFLS)
		electrostaticFactor = const_cast<LatticeAtomDSPER *>(this)->getGFLSFactor();
	if(ener < 0)
		ener = 0;
	return param->_sper[_type]._prefSPER[plane] * exp(-ener/kT) * electrostaticFactor;
}

/*** Obtains the plane when coordination is 1 */
unsigned LatticeAtomDSPER::getPlane(unsigned coord) const
{
	switch(coord)
	{
		case 4:
		case 3:
		case 2:
		{
			int coord2 = getCoordination(1) + coord;
			if(coord2 >= 10)
				return LatticeDiamondParam::P100_10;
			else if(coord2 == 9)
				return LatticeDiamondParam::P100_9;
			else if(coord2 == 8)
				return LatticeDiamondParam::P100_8;
			else if(coord2 == 7)
				return LatticeDiamondParam::P100_7;
			else if(coord2 <= 6)
				return LatticeDiamondParam::P100_6;
			break;
		}
		case 1:
		{
			//definition of 110: having a crystalline neighbor of a neighbor
			unsigned orient = 15;
			for(unsigned i=0; i<FIRSTN; ++i)
			{
				LatticeAtomDSPER *pNei = static_cast<LatticeAtomDSPER *>(_neighbors[i]);
				if(pNei && pNei->getPerformed()) //if "amorphous" neighbor
				{
					orient = pNei->_orientation;
					break;
				}
			}
			for(unsigned i=0; i<FIRSTN; ++i)
			{
				LatticeAtomDSPER *pNei = static_cast<LatticeAtomDSPER *>(_neighbors[i]);
				if(pNei && !pNei->getPerformed()) //if "amorphous" neighbor
				{
					unsigned howMany = 0;
					for(unsigned j=0; j<FIRSTN; ++j)
					{
						LatticeAtomDSPER *pNeiNei = static_cast<LatticeAtomDSPER *>(pNei->_neighbors[j]);
						if(pNeiNei && pNeiNei != this && pNeiNei->getPerformed() && pNeiNei->_orientation/2 == orient/2)
							howMany++;
					}
					if(howMany >= 1)
						return LatticeDiamondParam::P110;
					/*else if(howMany > 1)
						return LatticeDiamondParam::POTHER;*/
				}
			}
			// 111 or nothing
			for(unsigned i=FIRSTN; i<SECONDN; ++i)
			{
				LatticeAtomDSPER *pNei = static_cast<LatticeAtomDSPER *>(_neighbors[i]);
				if(pNei && !pNei->getPerformed()) //if "amorphous" neighbor
				{
					LatticeAtom *pSolid = 0;
					//unsigned coord=0;
					for(unsigned j=0; j<FIRSTN; ++j)
					{
						if(pNei->_neighbors[j] && pNei->_neighbors[j]->getPerformed())
					//	{
							pSolid = pNei->_neighbors[j];
							//coord++;
					//	}
					}
					if(/*coord == 1 && */pSolid && std::find(_neighbors+SECONDN, _neighbors+THIRDN, pSolid) != _neighbors+THIRDN)
						//be sure if forms an "hexagon"
						return LatticeDiamondParam::P111;
				}
			}
			return LatticeDiamondParam::PNONE;
		break;
		}
		default:
		break;
	}
	return LatticeDiamondParam::PNONE;
}


double LatticeAtomDSPER::getGFLSFactor()
{
	const LatticeDiamondParam  *param = static_cast<const LatticeDiamondParam *>(_pDomain->_pLaPar[_basicMat]);
	double                V   = 0.;
	double                Eg  = 0.;
	double                Eg0 = 0.;
	double                kT  = _pDomain->_pRM->getkT();

	const double          g_M  = param->_g_M;
	const double          g_0  = param->_g_0;
	const double          g_P  = param->_g_P;

	std::set<Electrostatics::MeshNode *> sMN;

	_pDomain->_pMesh->getParticleToNodeHandler()->getParticleNodes(static_cast<LatticeAtom *>(this), sMN);
	if (sMN.size() == 0) {
		return 1.;
	}

	for (std::set<Electrostatics::MeshNode *>::const_iterator it = sMN.begin(); it != sMN.end(); ++it) {
		V   += (*it)->getWeight(this) * (*it)->getPotential();
		Eg  += (*it)->getWeight(this) * (*it)->getBandgap();
		Eg0 += (*it)->getWeight(this) * (*it)->getBandgap0();
	}

	if (Eg <= 0. || Eg0 <= 0.) {
		return 1.;
	}

	const double Ef  = 0;
	const double VBM = - V - 0.5 * Eg;

    // Note: defect level scales with bandgap
	const double E_lkmc_M    = Ef - (VBM + Eg / Eg0 * param->_E_M0);           // Ef  - E(0,-)
	const double E_lkmc_Mi   = Ef - (- 0.5 * Eg + Eg / Eg0 * param->_E_M0);    // Efi - E(0,-)

	const double E_lkmc_P    = (VBM + Eg / Eg0 * param->_E_P0) - Ef;           // Ef  - E(+,0)
	const double E_lkmc_Pi   = (- 0.5 * Eg + Eg / Eg0 * param->_E_P0) - Ef;    // Efi - E(+,0)

	double M0d   = exp(E_lkmc_M / kT);
	double P0d   = exp(E_lkmc_P / kT);

	double M0i   = exp(E_lkmc_Mi / kT);
	double P0i   = exp(E_lkmc_Pi / kT);

	return (1. + g_M / g_0 * M0d + g_P / g_0 * P0d) / (1. + g_M / g_0 * M0i + g_P / g_0 * P0i);
}

unsigned LatticeAtomDSPER::idxMinVector(std::vector<double> &v) const
{
	if (v.size() <= 1)
		return 0;

	unsigned minIdx = 0;
	unsigned minVal = v[0];

	for (unsigned i = 1; i < v.size(); ++i) {
		if (v[i] < minVal) {
			minIdx = i;
			minVal = v[i];
		}
	}
	return minIdx;
}

unsigned LatticeAtomDSPER::idxMaxVector(std::vector<double> &v) const
{
	if (v.size() <= 1)
		return 0;

	unsigned maxIdx = 0;
	unsigned maxVal = v[0];

	for (unsigned i = 1; i < v.size(); ++i) {
		if (v[i] > maxVal) {
			maxIdx = i;
			maxVal = v[i];
		}
	}
	return maxIdx;
}

// calculate the cross product of v1 and v2
void LatticeAtomDSPER::cross_product(const ublas::vector<double> &v1, const ublas::vector<double> &v2, ublas::vector<double> &res) const{
	if (v1.size() != 3 && v2.size() != 3)
		return ;

	res = ublas::zero_vector<double>(v1.size());

	res(0) = v1(1) * v2(2) - v1(2) * v2(1);
	res(1) = v1(2) * v2(0) - v1(0) * v2(2);
	res(2) = v1(0) * v2(1) - v1(1) * v2(0);
}

void LatticeAtomDSPER::normalize(ublas::vector<double> &v) const
{
	float s = 0.;

	for (unsigned i = 0; i < v.size(); ++i) {
		if (fabs(v(i)) < 1.e-4)
			v(i) = 0.;
		else
			v(i) = fabs(v(i));
		s += fabs(v(i));
	}

	if (s != 0.)
		v /= s;
}

void LatticeAtomDSPER::direction110(ublas::vector<double> &d) const
{
	std::vector<ublas::vector<double> > c;

	ublas::vector<double> config110[2];
	ublas::vector<double> basis110[2];

	config110[0] = ublas::zero_vector<double>(3);
	config110[1] = ublas::zero_vector<double>(3);
	basis110[0]  = ublas::zero_vector<double>(3);
	basis110[1]  = ublas::zero_vector<double>(3);

	config110[0](0) = getCoordinates()._x;
	config110[0](1) = getCoordinates()._y;
	config110[0](2) = getCoordinates()._z;

	for (unsigned i = 0; i < FIRSTN; ++i)
	{
		LatticeAtomDSPER *pNei = static_cast<LatticeAtomDSPER *>(_neighbors[i]);

		if (pNei && pNei->getPerformed()) {
			basis110[0](0) = pNei->getCoordinates()._x;
			basis110[0](1) = pNei->getCoordinates()._y;
			basis110[0](2) = pNei->getCoordinates()._z;
		}

		if (pNei && !pNei->getPerformed()) //if "amorphous" neighbor
		{
			unsigned howMany = 0;

			for(unsigned j = 0; j < FIRSTN; ++j)
			{
				LatticeAtomDSPER *pNeiNei = static_cast<LatticeAtomDSPER *>(pNei->_neighbors[j]);
				if(pNeiNei && pNeiNei != this && pNeiNei->getPerformed()) {
					config110[1](0) = pNei->getCoordinates()._x;
					config110[1](1) = pNei->getCoordinates()._y;
					config110[1](2) = pNei->getCoordinates()._z;

					basis110[1](0)  = pNeiNei->getCoordinates()._x;
					basis110[1](1)  = pNeiNei->getCoordinates()._y;
					basis110[1](2)  = pNeiNei->getCoordinates()._z;

					howMany++;
				}
			}
		}
	}

	_pDomain->_pMesh->setPeriodicRelative(config110[0], config110[1]);

	ublas::vector<double> r = config110[0] + 0.5 * config110[1];

	_pDomain->_pMesh->setPeriodicRelative(r, basis110[0]);
	_pDomain->_pMesh->setPeriodicRelative(r, basis110[1]);

	d = basis110[0] + basis110[1];

	normalize(d);
}

void LatticeAtomDSPER::direction100(ublas::vector<double> &d) const
{
	std::vector<ublas::vector<double> > c; // crystalline 1st neighbors
	std::vector<double>                 n;

	ublas::vector<double> ref = ublas::zero_vector<double>(3);

	ref(0) = getCoordinates()._x;
	ref(1) = getCoordinates()._y;
	ref(2) = getCoordinates()._z;

	for(unsigned i = 0; i < FIRSTN; ++i)
		if(_neighbors[i]) {
			ublas::vector<double> v = ublas::zero_vector<double>(3);
			v(0) = _neighbors[i]->getCoordinates()._x;
			v(1) = _neighbors[i]->getCoordinates()._y;
			v(2) = _neighbors[i]->getCoordinates()._z;

			_pDomain->_pMesh->setPeriodicRelative(ref, v);

			if(_neighbors[i]->getPerformed())
			{
				c.push_back(v);
				unsigned k = 0;
				for (unsigned j = 0; j < SECONDN; ++j)
					if(static_cast<LatticeAtomDSPER *>(_neighbors[i])->_neighbors[j])
						if (static_cast<LatticeAtomDSPER *>(_neighbors[i])->_neighbors[j]->getPerformed())
							++k;
				n.push_back(k);
			}
		}

	unsigned coord = c.size();

	if (coord >= 2) {
		unsigned              i;
		ublas::vector<double> c1;
		ublas::vector<double> c2;
		ublas::vector<double> d100;

		i  = idxMaxVector(n);
		c1 = c[i];
		c.erase(c.begin() + i);
		n.erase(n.begin() + i);
		i  = idxMaxVector(n);
		c2 = c[i];

		d100 = c1 + c2;
		normalize(d100);

		d = d100;
	}
	normalize(d);
}

double LatticeAtomDSPER::stressCorrection_2ndOrder(unsigned plane, unsigned coord, float kT) const
{
	const LatticeDiamondParam *param = static_cast<const LatticeDiamondParam *>(_pDomain->_pLaPar[_basicMat]);

	ublas::matrix<double> stress = ublas::zero_matrix<double>(3, 3);

	const double omega   = Domains::global()->PM()->getMaterial(_pElement->getMaterial())._molarVolume * 1.e-6;   // m3/mol
	const double Na      = AVOGADRO;   // mol
	const double q       = Q;          // eV

	if (_pElement->stress().empty())
		return (0.);

	stress(0, 0) = _pElement->stress(0);
	stress(1, 1) = _pElement->stress(1);
	stress(2, 2) = _pElement->stress(2);

	if (plane >= LatticeDiamondParam::P100_6 && plane <= LatticeDiamondParam::P100_10) {
		ublas::vector<double> d = ublas::zero_vector<double>(3);

		if (plane >= LatticeDiamondParam::P100_6 && plane <= LatticeDiamondParam::P100_10)
			direction100(d);
		if (plane == LatticeDiamondParam::P110)
			direction110(d);


		ublas::matrix<double> dv        = ublas::zero_matrix<double>(3, 3);
		ublas::matrix<double> kronecker = ublas::identity_matrix<double>(3, 3);

		double dvpar  = 0.;
		double dvperp = 0.;

		if (coord == 2) {
			dvpar  = param->_dvpar100_2;
			dvperp = param->_dvperp100_2;
		}
		else if (coord == 3) {
			dvpar  = param->_dvpar100_3;
			dvperp = param->_dvperp100_3;
		}

		if (plane == LatticeDiamondParam::P110) {
			dvpar  = param->_dvpar110;
			dvperp = param->_dvperp110;
		}

		for (unsigned i = 0; i < 3; ++i) {
			for (unsigned j = 0; j < 3; ++j) {
				dv(i, j) = dvpar * d(i) * d(j) + dvperp * (kronecker(i, j) - d(i) * d(j));
			}
		}

		ublas::matrix<double> tmp = ublas::prod(dv, stress);
		double E = 0.;

		for(unsigned k = 0; k < tmp.size1(); ++k)
			for(unsigned l = 0; l < tmp.size2(); ++l)
				E += tmp(k, l);
		E *= omega / Na / q;
		return E;
	}

	return 0.;
}

} /* namespace OKMC */
