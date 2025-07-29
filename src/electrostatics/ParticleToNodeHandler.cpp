/*
 * Author: Benoit Sklenard benoit.sklenard@cea.fr 
 * 
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

#include <cmath>
#include <vector>

#include "io/Diagnostic.h"
#include "io/ParameterManager.h"

#include "okmc/Defect.h"
#include "kernel/ParticleType.h"
#include "okmc/MobileParticleParam.h"

#include <boost/numeric/ublas/io.hpp>

#include "kernel/Domain.h"
#include "kernel/Mesh.h"

#include "ParticleToNodeHandler.h"

#include <ctime>

using namespace Electrostatics;
using namespace Kernel;
using namespace boost::numeric;

ParticleToNodeHandler::ParticleToNodeHandler(Kernel::Domain *pDomain, Kernel::Mesh *pMesh) {
	LOWMSG("Loading Particle To Node handler");
	_pMesh   = pMesh;
	_pDomain = pDomain;
}

ParticleToNodeHandler::~ParticleToNodeHandler() {

}

inline double ParticleToNodeHandler::getOverlap(double u1, double u2) {
	return 15./16. * ((u2 - u1) - 2./3. * (u2*u2*u2 - u1*u1*u1) + 1./5. * (u2*u2*u2*u2*u2 - u1*u1*u1*u1*u1));
}

void ParticleToNodeHandler::remove(OKMC::Particle *pPart) {
	std::map<OKMC::Particle *, std::set<MeshNode *> >::iterator mit;

	mit = _syncOKMC.find(pPart);
	if (mit != _syncOKMC.end()) {
		for (std::set<MeshNode *>::iterator sit = mit->second.begin(); sit != mit->second.end(); ++sit)
			(*sit)->remove(pPart);
		_syncOKMC.erase(mit);
	}
}

void ParticleToNodeHandler::remove(LKMC::LatticeAtom *pLA) {
	std::map<LKMC::LatticeAtom *, std::set<MeshNode *> >::iterator mit;

	mit = _syncLKMC.find(pLA);
	if (mit != _syncLKMC.end()) {
		for (std::set<MeshNode *>::iterator sit = mit->second.begin(); sit != mit->second.end(); ++sit)
			(*sit)->remove(pLA);
		_syncLKMC.erase(mit);
	}
}

double ParticleToNodeHandler::addWeightNode(MeshNode *pNode, LKMC::LatticeAtom *pLA) {
	M_TYPE mt = pLA->getElement()->getMaterial();
	double dr = _pDomain->_pLaPar[mt]->_orbitalRadius;

	ublas::vector<double> n  = ublas::zero_vector<double>(3);
	ublas::vector<double> c1 = ublas::zero_vector<double>(3);
	ublas::vector<double> c2 = ublas::zero_vector<double>(3);
	MeshNode ***pNodes = _pMesh->getNodes();

	n(0)  = pLA->getCoordinates()._x;
	n(1)  = pLA->getCoordinates()._y;
	n(2)  = pLA->getCoordinates()._z;

	c1(0) = pNode->_xm;
	c1(1) = pNode->_ym;
	c1(2) = pNode->_zm;
	c2(0) = pNode->_xp;
	c2(1) = pNode->_yp;
	c2(2) = pNode->_zp;

	_pMesh->setPeriodicRelative(n, c1); // c1 = c1 - n
	_pMesh->setPeriodicRelative(n, c2); // c2 = c2 - n

	c1 /= dr;
	c2 /= dr;

	if (c1(0) <= -1) c1(0) = -1.;
	if (c2(0) <= -1) c2(0) = -1.;
	if (c1(0) >= 1)  c1(0) = 1.;
	if (c2(0) >= 1)  c2(0) = 1.;

	if (c1(1) <= -1) c1(1) = -1.;
	if (c2(1) <= -1) c2(1) = -1.;
	if (c1(1) >= 1)  c1(1) = 1.;
	if (c2(1) >= 1)  c2(1) = 1.;

	if (c1(2) <= -1) c1(2) = -1.;
	if (c2(2) <= -1) c2(2) = -1.;
	if (c1(2) >= 1)  c1(2) = 1.;
	if (c2(2) >= 1)  c2(2) = 1.;

	const double wx  = getOverlap(c1(0), c2(0));
	const double wy  = getOverlap(c1(1), c2(1));
	const double wz  = getOverlap(c1(2), c2(2));

	const double w   = wx * wy * wz;

	if (w > 0) {
		unsigned i = pNode->_ix;
		unsigned j = pNode->_iy;
		unsigned k = pNode->_iz;

		if (_pMesh->getPeriodicX() && pNode->_ix == (_pMesh->getnx() - 1))
			i = 0;
		if (_pMesh->getPeriodicY() && pNode->_iy == (_pMesh->getny() - 1))
			j = 0;
		if (_pMesh->getPeriodicZ() && pNode->_iz == (_pMesh->getnz() - 1))
			k = 0;

		pNodes[i][j][k].insert(pLA, w);

	    std::pair<std::map<LKMC::LatticeAtom *, std::set<MeshNode *> >::iterator, bool > r;
	    r = _syncLKMC.insert(std::pair<LKMC::LatticeAtom *, std::set<MeshNode *> >(pLA, std::set<MeshNode *>()));
	    r.first->second.insert(&pNodes[i][j][k]);
	}

	return w;
}

double ParticleToNodeHandler::addWeightNode(MeshNode *pNode, OKMC::Particle *pPart) {
	M_TYPE mt = pPart->getElement()->getMaterial();
	P_TYPE pt = pPart->getPType();

	double dr = _pDomain->_pMPPar->_orbitalRadius[mt][pt];
	// double dr = 1.; // _pDomain->_pLaPar[mt]->_orbitalRadius;

	ublas::vector<double> n  = ublas::zero_vector<double>(3);
	ublas::vector<double> c1 = ublas::zero_vector<double>(3);
	ublas::vector<double> c2 = ublas::zero_vector<double>(3);
	MeshNode ***pNodes = _pMesh->getNodes();

	n(0)  = pPart->getCoordinates()._x;
	n(1)  = pPart->getCoordinates()._y;
	n(2)  = pPart->getCoordinates()._z;

	c1(0) = pNode->_xm;
	c1(1) = pNode->_ym;
	c1(2) = pNode->_zm;
	c2(0) = pNode->_xp;
	c2(1) = pNode->_yp;
	c2(2) = pNode->_zp;

	_pMesh->setPeriodicRelative(n, c1); // c1 = c1 - n
	_pMesh->setPeriodicRelative(n, c2); // c2 = c2 - n

	c1 /= dr;
	c2 /= dr;

	if (c1(0) <= -1) c1(0) = -1.;
	if (c2(0) <= -1) c2(0) = -1.;
	if (c1(0) >= 1)  c1(0) = 1.;
	if (c2(0) >= 1)  c2(0) = 1.;

	if (c1(1) <= -1) c1(1) = -1.;
	if (c2(1) <= -1) c2(1) = -1.;
	if (c1(1) >= 1)  c1(1) = 1.;
	if (c2(1) >= 1)  c2(1) = 1.;

	if (c1(2) <= -1) c1(2) = -1.;
	if (c2(2) <= -1) c2(2) = -1.;
	if (c1(2) >= 1)  c1(2) = 1.;
	if (c2(2) >= 1)  c2(2) = 1.;

	const double wx  = getOverlap(c1(0), c2(0));
	const double wy  = getOverlap(c1(1), c2(1));
	const double wz  = getOverlap(c1(2), c2(2));

	const double w   = wx * wy * wz;

	if (w > 0) {
		unsigned i = pNode->_ix;
		unsigned j = pNode->_iy;
		unsigned k = pNode->_iz;

		if (_pMesh->getPeriodicX() && pNode->_ix == (_pMesh->getnx() - 1))
			i = 0;
		if (_pMesh->getPeriodicY() && pNode->_iy == (_pMesh->getny() - 1))
			j = 0;
		if (_pMesh->getPeriodicZ() && pNode->_iz == (_pMesh->getnz() - 1))
			k = 0;

		pNodes[i][j][k].insert(pPart, w);

	    std::pair<std::map<OKMC::Particle *, std::set<MeshNode *> >::iterator, bool > r;
	    r = _syncOKMC.insert(std::pair<OKMC::Particle *, std::set<MeshNode *> >(pPart, std::set<MeshNode *>()));
	    r.first->second.insert(&pNodes[i][j][k]);
	}

	return w;
}

void ParticleToNodeHandler::insert(LKMC::LatticeAtom *pLA) {
	std::set<MeshNode*> MN;
	MeshNode* pMN = NULL;

	M_TYPE mt = pLA->getElement()->getMaterial();

	double dr = _pDomain->_pLaPar[mt]->_orbitalRadius;

	remove(pLA);
	MeshNode ***pNodes  = _pMesh->getNodes();
	pMN = _pMesh->getFirstNodeFromElement(pLA->getElement());


	const size_t dx = ceil(dr / fabs(pNodes[0][0][0]._xm - pNodes[0][0][0]._xp)) + 1;
	const size_t dy = ceil(dr / fabs(pNodes[0][0][0]._ym - pNodes[0][0][0]._yp)) + 1;
	const size_t dz = ceil(dr / fabs(pNodes[0][0][0]._zm - pNodes[0][0][0]._zp)) + 1;

	const size_t startx = (dx <= pMN->_ix ? pMN->_ix - dx : _pMesh->getnx() - 1 - dx + pMN->_ix);
	const size_t endx   = (pMN->_ix + dx) % _pMesh->getnx();
	const size_t starty = (dy <= pMN->_iy ? pMN->_iy - dy : _pMesh->getny() - 1 - dy + pMN->_iy);
	const size_t endy   = (pMN->_iy + dy) % _pMesh->getny();
	const size_t startz = (dz <= pMN->_iz ? pMN->_iz - dz : _pMesh->getnz() - 1 - dz + pMN->_iz);
	const size_t endz   = (pMN->_iz + dz) % _pMesh->getnz();

	double w = 0.;
	for (size_t ix = startx; ix != endx; ix = (ix + 1) % _pMesh->getnx()) {
		for (size_t iy = starty; iy != endy; iy = (iy + 1) % _pMesh->getny()) {
			for (size_t iz = startz; iz != endz; iz = (iz + 1) % _pMesh->getnz()) {
				w += addWeightNode(&pNodes[ix][iy][iz], pLA);
			}
		}
	}

	if (_syncLKMC.find(pLA) ==_syncLKMC.end())
		ERRORMSG("LatticeAtom has not been inserted!");
}


void ParticleToNodeHandler::insert(OKMC::Particle *pPart) {
	std::set<MeshNode*> MN;
	MeshNode* pMN = NULL;

	M_TYPE mt = pPart->getElement()->getMaterial();
	P_TYPE pt = pPart->getPType();

	if (!_pDomain->_pMPPar->_mapToGrid[mt][pt]) {
	//	WARNINGMSG("Not inserted (OKMC) " << Domains::global()->PM()->getMaterialName(mt) << " " << Domains::global()->PM()->getParticleName(mt, pt));
		return ;
	}

	double dr = _pDomain->_pMPPar->_orbitalRadius[mt][pt];

	remove(pPart);
	MeshNode ***pNodes  = _pMesh->getNodes();
	pMN = _pMesh->getFirstNodeFromElement(pPart->getElement());

	const size_t dx = ceil(dr / fabs(pNodes[0][0][0]._xm - pNodes[0][0][0]._xp)) + 1;
	const size_t dy = ceil(dr / fabs(pNodes[0][0][0]._ym - pNodes[0][0][0]._yp)) + 1;
	const size_t dz = ceil(dr / fabs(pNodes[0][0][0]._zm - pNodes[0][0][0]._zp)) + 1;

	const size_t startx = (dx <= pMN->_ix ? pMN->_ix - dx : _pMesh->getnx() - 1 - dx + pMN->_ix);
	const size_t endx   = (pMN->_ix + dx) % _pMesh->getnx();
	const size_t starty = (dy <= pMN->_iy ? pMN->_iy - dy : _pMesh->getny() - 1 - dy + pMN->_iy);
	const size_t endy   = (pMN->_iy + dy) % _pMesh->getny();
	const size_t startz = (dz <= pMN->_iz ? pMN->_iz - dz : _pMesh->getnz() - 1 - dz + pMN->_iz);
	const size_t endz   = (pMN->_iz + dz) % _pMesh->getnz();

	double w = 0.;
	for (size_t ix = startx; ix != endx; ix = (ix + 1) % _pMesh->getnx()) {
		for (size_t iy = starty; iy != endy; iy = (iy + 1) % _pMesh->getny()) {
			for (size_t iz = startz; iz != endz; iz = (iz + 1) % _pMesh->getnz()) {
				w += addWeightNode(&pNodes[ix][iy][iz], pPart);
			}
		}
	}

	if (_syncOKMC.find(pPart) ==_syncOKMC.end())
		ERRORMSG("Particle has not been inserted!");
}

void ParticleToNodeHandler::getParticleNodes(OKMC::Particle *pPart, std::set<MeshNode *> &nodes) {
	std::map<OKMC::Particle *, std::set<MeshNode *> >::iterator mit = _syncOKMC.find(pPart);

	if (mit != _syncOKMC.end())
		for (std::set<MeshNode *>::iterator sit = mit->second.begin(); sit != mit->second.end(); ++sit)
			nodes.insert(*sit);
}

void ParticleToNodeHandler::getParticleNodes(LKMC::LatticeAtom *pLA, std::set<MeshNode *> &nodes) {
	std::map<LKMC::LatticeAtom *, std::set<MeshNode *> >::iterator mit = _syncLKMC.find(pLA);

	if (mit != _syncLKMC.end())
		for (std::set<MeshNode *>::iterator sit = mit->second.begin(); sit != mit->second.end(); ++sit)
			nodes.insert(*sit);
}

void ParticleToNodeHandler::print() {
	LOWMSG(_syncOKMC.size() << " particles in ParticleToNodeHandler:");

	if (!_syncOKMC.empty()) {
		unsigned i = 0;

		for (std::map<OKMC::Particle *, std::set<MeshNode *> >::iterator mit = _syncOKMC.begin(); mit != _syncOKMC.end(); ++mit) {
			double w = 0.;
			for (std::set<MeshNode *>::iterator sit = mit->second.begin(); sit != mit->second.end(); ++sit) {
				std::map<OKMC::Particle *, double>::iterator mit2 = (*sit)->_mPart.find(mit->first);
				if (mit2 != (*sit)->_mPart.end())
					w += mit2->second;
				else
					LOWMSG("Got a null value!");
			}
			LOWMSG("Particle " << i++ << " " <<  mit->first->getCoordinates() << " is affected to " << mit->second.size() << " nodes" << " weight= " << w);
		}
	}
}

void ParticleToNodeHandler::printLKMC() {
	LOWMSG(_syncLKMC.size() << " LKMC::LatticeAtom in ParticleToNodeHandler:");

	if (!_syncLKMC.empty()) {
		unsigned i = 0;

		for (std::map<LKMC::LatticeAtom *, std::set<MeshNode *> >::iterator mit = _syncLKMC.begin(); mit != _syncLKMC.end(); ++mit) {
			double w = 0.;
			for (std::set<MeshNode *>::iterator sit = mit->second.begin(); sit != mit->second.end(); ++sit) {
				std::map<LKMC::LatticeAtom *, double>::iterator mit2 = (*sit)->_mLA.find(mit->first);
				if (mit2 != (*sit)->_mLA.end())
					w += mit2->second;
				else
					LOWMSG("Got a null value!");
			}
				LOWMSG("LKMC::LaticeAtom " << ++i << " " <<  mit->first->getCoordinates() << " is affected to " << mit->second.size() << " nodes" << " weight= " << w);
		}
	}
}
