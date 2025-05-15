/*
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

#include "Mesh.h"
#include "lkmc/LatticeSite.h"
#include "lkmc/LatticeAtom.h"
#include "domains/MCClient.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "Domain.h"
#include "SubDomain.h"
#include "MeshParam.h"
#include "RateManager.h"
#include "okmc/Particle.h"
#include "okmc/Interface.h"
#include "okmc/InterfaceParam.h"
#include "okmc/MobileParticle.h"
#include "okmc/MobileParticleParam.h"
#include "domains/Splitter.h"
#include "Constants.h"
#include <cmath>
#include <fstream>

using std::vector;
using std::pair;

using namespace Electrostatics;
using namespace Kernel;
using LKMC::LatticeSite;

Mesh::Mesh(Domain *p, const Coordinates &m, const Coordinates &M,
const vector<float> &x, const vector<float> &y, const vector<float> &z,
const Domains::MCClient *pCli) : _min(m), _max(M)
{
	_pDomain = p;
	_xlines = x;
	_ylines = y;
	_zlines = z;
	_zCells  = _zlines.size() -1;
	_yzCells = (_zlines.size() -1)*(_ylines.size() -1);
	unsigned NElements = (_xlines.size() -1) * (_ylines.size() -1) * (_zlines.size() -1);
	_elements.assign(NElements,MeshElement(_pDomain));
	unsigned idx=0;
	for(vector<MeshElement>::iterator it=_elements.begin(); it!=_elements.end(); ++it)
	{                            
		it->_index = idx++;
		Coordinates center;
		getCenter(center, idx-1);
		it->_mat = pCli->getMaterial()->operator ()(center, it->_index);
		Kernel::Coordinates m, M;
		getCorners(idx-1, m, M);
		it->_subDomainIdx = Domains::global()->getSplitter()->getSubDomain(m, M);
		it->_subDomainLevel = Domains::global()->getSplitter()->getLevel(m, M);
		it->_amorphParts = 0.;
		float vol = (M._x - m._x)*(M._y - m._y)*(M._z - m._z);
		if(_pDomain->_pAlPar->_isAlloy[it->_mat])
			it->_AAtoms = round(Domains::global()->PM()->getMaterial(it->_mat)._densityAlloyCm3 * vol * 1e-21);
		else
			it->_AAtoms = round(Domains::global()->PM()->getMaterial(it->_mat)._densityCm3 * vol * 1e-21);
	}
	for(unsigned i=0; i< _ylines.size() -1; ++i) _yParticles.push_back(0);
	for(unsigned i=0; i< _xlines.size() -1; ++i) _xParticles.push_back(0);

	_periodicX = _pDomain->_pMePar->_periodicX;
	_periodicY = _pDomain->_pMePar->_periodicY;
	_periodicZ = _pDomain->_pMePar->_periodicZ;
	_xsize = _xlines.back() - _xlines.front();
	_ysize = _ylines.back() - _ylines.front();
	_zsize = _zlines.back() - _zlines.front();

	_poissonX  = _pDomain->_pMePar->_poissonX;
	_poissonY  = _pDomain->_pMePar->_poissonY;
	_poissonZ  = _pDomain->_pMePar->_poissonZ;

	_potentialX = _pDomain->_pMePar->_potentialX;
	_potentialY = _pDomain->_pMePar->_potentialY;
	_potentialZ = _pDomain->_pMePar->_potentialZ;

	_pPNH = new ParticleToNodeHandler(_pDomain, this);

	_pNodes = NULL;

	for(M_TYPE mt=0; mt<Domains::global()->PM()->getNMaterials(); ++mt) {
		for(uint32_t x = 1u; x < _xlines.size(); ++x) {
            for(uint32_t y = 1u; y < _ylines.size(); ++y) {
				auto const xm1 = x - 1u;
				auto const ym1 = y - 1u;
				_longHopFactor[mt][xm1 * (_ylines.size() - 1u) + ym1] =
std::max(unsigned(1), unsigned(std::min(_xlines[x] - _xlines[xm1], _ylines[y] - _ylines[ym1]) / _pDomain->_pMePar->_lambda[mt]));
			}
		}
	}
	buildNodes();
	initNodes();
	buildNeighbors();

	print();
	check();
}

Mesh::~Mesh()
{
	for (unsigned ix = 0; ix < _xlines.size(); ++ix) {
		for (unsigned iy = 0; iy < _ylines.size(); ++iy)
			delete [] _pNodes[ix][iy];
		delete [] _pNodes[ix];
	}
	delete [] _pNodes;

	delete _pPNH;
}

void Mesh::buildNodes() {
	unsigned idxActive   = 0;
	unsigned idxInactive = -1;

	LOWMSG2("Building nodes... ");

	_pNodes = new MeshNode **[_xlines.size()];
	for (unsigned ix = 0; ix < _xlines.size(); ++ix) {
		_pNodes[ix] = new MeshNode *[_ylines.size()];
		for (unsigned iy = 0; iy < _ylines.size(); ++iy) {
			_pNodes[ix][iy] = new MeshNode [_zlines.size()];
			for (unsigned iz = 0; iz < _zlines.size(); ++iz) {
				_pNodes[ix][iy][iz]._active = true;
				_pNodes[ix][iy][iz]._ix     = ix;
				_pNodes[ix][iy][iz]._iy     = iy;
				_pNodes[ix][iy][iz]._iz     = iz;

				if (_periodicX && ix == (_xlines.size() - 1))
					_pNodes[ix][iy][iz]._active = false;
				if (_periodicY && iy == (_ylines.size() - 1))
					_pNodes[ix][iy][iz]._active = false;
				if (_periodicZ && iz == (_zlines.size() - 1))
					_pNodes[ix][iy][iz]._active = false;

				if (_poissonX == "Dirichlet" && (ix == 0 || ix == (_xlines.size() - 1))) {
					_pNodes[ix][iy][iz]._active    = false;
					_pNodes[ix][iy][iz]._V = _potentialX;
				}
				if (_poissonY == "Dirichlet" && (iy == 0 || iy == (_ylines.size() - 1))) {
					_pNodes[ix][iy][iz]._active    = false;
					_pNodes[ix][iy][iz]._V = _potentialY;
				}
				if (_poissonZ == "Dirichlet" && (iz == 0 || iz == (_zlines.size() - 1))) {
					_pNodes[ix][iy][iz]._active    = false;
					_pNodes[ix][iy][iz]._V = _potentialZ;
				}

				if (_pNodes[ix][iy][iz]._active)
					_pNodes[ix][iy][iz]._index   = idxActive++;
				else
					_pNodes[ix][iy][iz]._index   = idxInactive--;

				_pNodes[ix][iy][iz]._coord   = Kernel::Coordinates(_xlines[ix], _ylines[iy], _zlines[iz]);
			}
		}
	}

	LOWMSG("Done");
}

void Mesh::initNodes() {
	for (unsigned ix = 0; ix < _xlines.size(); ++ix) {
		for (unsigned iy = 0; iy < _ylines.size(); ++iy) {
			for (unsigned iz = 0; iz < _zlines.size(); ++iz) {
				getMidX(ix, iy, iz, _pNodes[ix][iy][iz]._xm, _pNodes[ix][iy][iz]._xp);
				getMidY(ix, iy, iz, _pNodes[ix][iy][iz]._ym, _pNodes[ix][iy][iz]._yp);
				getMidZ(ix, iy, iz, _pNodes[ix][iy][iz]._zm, _pNodes[ix][iy][iz]._zp);

				if(_pNodes[ix][iy][iz]._active) {
					_pNodes[ix][iy][iz]._volume += (_pNodes[ix][iy][iz]._xp - _pNodes[ix][iy][iz]._xm) *
	                                               (_pNodes[ix][iy][iz]._yp - _pNodes[ix][iy][iz]._ym) *
	                                               (_pNodes[ix][iy][iz]._zp - _pNodes[ix][iy][iz]._zm);
				}
				else {
					unsigned i = ix;
					unsigned j = iy;
					unsigned k = iz;
					if (_periodicX && ix == (_xlines.size() - 1))
						i = 0;
					if (_periodicY && iy == (_ylines.size() - 1))
						j = 0;
					if (_periodicZ && iz == (_zlines.size() - 1))
						k = 0;

					_pNodes[i][j][k]._volume += (_pNodes[ix][iy][iz]._xp - _pNodes[ix][iy][iz]._xm) *
                            					(_pNodes[ix][iy][iz]._yp - _pNodes[ix][iy][iz]._ym) *
                                                (_pNodes[ix][iy][iz]._zp - _pNodes[ix][iy][iz]._zm);
				}
			}
		}
	}

	for (unsigned ix = 0; ix < _xlines.size(); ++ix) {
		for (unsigned iy = 0; iy < _ylines.size(); ++iy) {
			for (unsigned iz = 0; iz < _zlines.size(); ++iz) {
				if(_pNodes[ix][iy][iz]._active) {
					double xm = _pNodes[ix][iy][iz]._xm;
					double xp = _pNodes[ix][iy][iz]._xp;
					double ym = _pNodes[ix][iy][iz]._ym;
					double yp = _pNodes[ix][iy][iz]._yp;
					double zm = _pNodes[ix][iy][iz]._zm;
					double zp = _pNodes[ix][iy][iz]._zp;

					if (_periodicX && ix == 0)
						xm = _pNodes[_xlines.size() - 1][iy][iz]._xm;
					if (_periodicY && iy == 0)
						ym = _pNodes[ix][_ylines.size() - 1][iz]._ym;
					if (_periodicZ && iz == 0)
						zm = _pNodes[ix][iy][_zlines.size() - 1]._zm;

					MeshElement *pME[8];

					unsigned xmin = 0;
					unsigned xmax = _xlines.size() - 1;
					unsigned ymin = 0;
					unsigned ymax = _ylines.size() - 1;
					unsigned zmin = 0;
					unsigned zmax = _zlines.size() - 1;

					for (unsigned i = 0; i < 8; ++i)
						pME[i] = NULL;

					if ((_periodicX || ix != xmin) && (_periodicY || iy != ymin) && (_periodicZ || iz != zmin))
						pME[0] = getElement(getIndexFromCoordinates(Coordinates(xm, ym, zm)));
					if ((_periodicX || ix != xmax) && (_periodicY || iy != ymin) && (_periodicZ || iz != zmin))
						pME[1] = getElement(getIndexFromCoordinates(Coordinates(xp, ym, zm)));
					if ((_periodicX || ix != xmin) && (_periodicY || iy != ymax) && (_periodicZ || iz != zmin))
						pME[2] = getElement(getIndexFromCoordinates(Coordinates(xm, yp, zm)));
					if ((_periodicX || ix != xmin) && (_periodicY || iy != ymin) && (_periodicZ || iz != zmax))
						pME[3] = getElement(getIndexFromCoordinates(Coordinates(xm, ym, zp)));
					if ((_periodicX || ix != xmin) && (_periodicY || iy != ymax) && (_periodicZ || iz != zmax))
						pME[4] = getElement(getIndexFromCoordinates(Coordinates(xm, yp, zp)));
					if ((_periodicX || ix != xmax) && (_periodicY || iy != ymin) && (_periodicZ || iz != zmax))
						pME[5] = getElement(getIndexFromCoordinates(Coordinates(xp, ym, zp)));
					if ((_periodicX || ix != xmax) && (_periodicY || iy != ymax) && (_periodicZ || iz != zmin))
						pME[6] = getElement(getIndexFromCoordinates(Coordinates(xp, yp, zm)));
					if ((_periodicX || ix != xmax) && (_periodicY || iy != ymax) && (_periodicZ || iz != zmax))
						pME[7] = getElement(getIndexFromCoordinates(Coordinates(xp, yp, zp)));

				    std::pair<std::map<MeshNode *, std::set<MeshElement *> >::iterator, bool > r1;
				    std::pair<std::map<MeshElement *, std::set<MeshNode *> >::iterator, bool > r2;

				    r1 = _MN2ME.insert(std::pair<MeshNode *, std::set<MeshElement *> >(&_pNodes[ix][iy][iz], std::set<MeshElement *>()));
					for (unsigned i = 0; i < 8; ++i) {
						if (pME[i]) {
							r1.first->second.insert(pME[i]);
						    r2 = _ME2MN.insert(std::pair<MeshElement *, std::set<MeshNode *> >(pME[i], std::set<MeshNode *>()));
						    r2.first->second.insert(&_pNodes[ix][iy][iz]);
						}
					}
				}
			}
		}
	}
}

void Mesh::getMidX(unsigned ix, unsigned iy, unsigned iz, double &xm, double &xp) {
	if (ix != 0 && ix != (_xlines.size() - 1)) {
		xm = 1. / 2. * (_pNodes[ix][iy][iz]._coord._x + _pNodes[ix-1][iy][iz]._coord._x);
		xp = 1. / 2. * (_pNodes[ix][iy][iz]._coord._x + _pNodes[ix+1][iy][iz]._coord._x);
	} else {
		if (ix == 0) {
			xm = _pNodes[ix][iy][iz]._coord._x;
			xp = 1. / 2. * (_pNodes[ix][iy][iz]._coord._x + _pNodes[ix+1][iy][iz]._coord._x);
		}
		if (ix == (_xlines.size() - 1)) {
			xm = 1. / 2. * (_pNodes[ix][iy][iz]._coord._x + _pNodes[ix-1][iy][iz]._coord._x);
			xp = _pNodes[ix][iy][iz]._coord._x;
		}
	}
}

void Mesh::getMidY(unsigned ix, unsigned iy, unsigned iz, double &ym, double &yp) {
	if (iy != 0 && iy != (_ylines.size() - 1)) {
		ym = 1. / 2. * (_pNodes[ix][iy][iz]._coord._y + _pNodes[ix][iy-1][iz]._coord._y);
		yp = 1. / 2. * (_pNodes[ix][iy][iz]._coord._y + _pNodes[ix][iy+1][iz]._coord._y);
	} else {
		if (iy == 0) {
			ym = _pNodes[ix][iy][iz]._coord._y;
			yp = 1. / 2. * (_pNodes[ix][iy][iz]._coord._y + _pNodes[ix][iy+1][iz]._coord._y);
		}
		if (iy == (_ylines.size() - 1)) {
			ym = 1. / 2. * (_pNodes[ix][iy][iz]._coord._y + _pNodes[ix][iy-1][iz]._coord._y);
			yp = _pNodes[ix][iy][iz]._coord._y;
		}
	}
}

void Mesh::getMidZ(unsigned ix, unsigned iy, unsigned iz, double &zm, double &zp) {
	if (iz != 0 && iz != (_zlines.size() - 1)) {
		zm = 1 / 2. * (_pNodes[ix][iy][iz]._coord._z + _pNodes[ix][iy][iz-1]._coord._z);
		zp = 1 / 2. * (_pNodes[ix][iy][iz]._coord._z + _pNodes[ix][iy][iz+1]._coord._z);
	} else {
		if (iz == 0) {
			zm = _pNodes[ix][iy][iz]._coord._z;
			zp = 1. / 2. * (_pNodes[ix][iy][iz]._coord._z + _pNodes[ix][iy][iz+1]._coord._z);
		}
		if (iz == (_zlines.size() - 1)) {
			zm = 1. / 2. * (_pNodes[ix][iy][iz]._coord._z + _pNodes[ix][iy][iz-1]._coord._z);
			zp = _pNodes[ix][iy][iz]._coord._z;
		}
	}
}

void Mesh::getElementsFromNode(MeshNode *pMN, std::set<MeshElement *> &ME) const {
	std::map<MeshNode *, std::set<MeshElement *> >::const_iterator it = _MN2ME.find(pMN);

	if (it != _MN2ME.end()) {
		ME = it->second;
	}
}

void Mesh::getNodesFromElement(MeshElement *pME, std::set<MeshNode*> &MN) const {
	std::map<MeshElement *, std::set<MeshNode *> >::const_iterator it = _ME2MN.find(pME);

	if (it != _ME2MN.end()) {
		MN = it->second;
	}
}

MeshNode * Mesh::getFirstNodeFromElement(MeshElement *pME) const {
	std::map<MeshElement *, std::set<MeshNode *> >::const_iterator it = _ME2MN.find(pME);
	return *(it->second.begin());
}

/*
 * FIXME: isDirichletX, isDirichletY and isDirichletZ will be replaced *soon*
 * std::string comparison should not be used here but it was fast to implement! :-)
 */
bool Mesh::isDirichletX() const {
	if (_poissonX == "Dirichlet")
		return true;
	else
		return false;
}

bool Mesh::isDirichletY() const {
	if (_poissonY == "Dirichlet")
		return true;
	else
		return false;
}

bool Mesh::isDirichletZ() const {
	if (_poissonZ == "Dirichlet")
		return true;
	else
		return false;
}

void Mesh::updatePotential(std::set<MeshNode *> &nodes) {
	std::set<MeshElement *> elements;

	for (std::set<MeshNode *>::iterator sit = nodes.begin(); sit != nodes.end(); ++sit) {
		std::map<MeshNode *, std::set<MeshElement *> >::iterator mit = _MN2ME.find(*sit);

		for (std::set<MeshElement *>::iterator it = mit->second.begin(); it != mit->second.end(); ++it)
			elements.insert(*it);
	}

	for (std::set<MeshElement *>::iterator sit = elements.begin(); sit != elements.end(); ++sit) {
		std::map<MeshElement *, std::set<MeshNode *> >::iterator mit = _ME2MN.find(*sit);

		if (mit != _ME2MN.end()) {
			mit->first->_V = 0;

			for (std::set<MeshNode *>::iterator it = mit->second.begin(); it != mit->second.end(); ++it) {
				mit->first->_V += (*it)->_V;
			}
			mit->first->_V /= static_cast<double>(mit->second.size());
		}
	}
}

void Mesh::updatePotential() {
	for(long unsigned i = 0; i < _elements.size(); ++i)
	{
		std::map<MeshElement *, std::set<MeshNode *> >::iterator mit = _ME2MN.find(&_elements[i]);

		if (mit != _ME2MN.end()) {
			mit->first->_V = 0;

			for (std::set<MeshNode *>::iterator sit = mit->second.begin(); sit != mit->second.end(); ++sit) {
				mit->first->_V += (*sit)->_V;
			}
			mit->first->_V /= static_cast<double>(mit->second.size());
		}
	}
}

void Mesh::print() const
{
	std::string oldLine;
	for(Kernel::M_TYPE mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
		LOWMSG(int(mt) << " -> " << Domains::global()->PM()->getMaterialName(mt));
	for(unsigned ix=0; ix<_xlines.size()-1; ++ix)
	{
		std::stringstream ss;
		for(unsigned iy=0; iy<_ylines.size()-1;++iy)
		{
			const MeshElement &me = _elements[getIndexFromIndices(ix, iy, 0)];
			if(me.getMaterial() < Kernel::MAX_MATERIALS)
				ss << int(me.getMaterial());
			else
				ss << "?";
		}
		if(oldLine != ss.str() || ix == _xlines.size() - 2)
		{
			LOWMSG2(ss.str());
			LOWMSG2(" - ");
			for(unsigned iy=0; iy<_ylines.size()-1;++iy)
			{
				const MeshElement &me = _elements[getIndexFromIndices(ix, iy, 0)];
				LKMC::LatticeAtom *pLA = dynamic_cast<LKMC::LatticeAtom *>(me.getFirstLS());
				LOWMSG2((pLA? 'X' : ' '));
			}

			LOWMSG2("  " << _xlines[ix] << ":" << _xlines[ix+1] << "\n");
		}
		oldLine = ss.str();
	}
}

void Mesh::printDomains() const
{
	LOWMSG("----------- SubDomains -- Levels --");
	for(unsigned ix=0; ix<_xlines.size()-1; ++ix)
	{
		for(unsigned iy=0; iy<_ylines.size()-1;++iy)
		{
			const MeshElement &me = _elements[getIndexFromIndices(ix, iy, 0)];
			LOWMSG2( char('A' + me._subDomainIdx));
		}
		LOWMSG2(" - ");
		for(unsigned iy=0; iy<_ylines.size()-1;++iy)
		{
			const MeshElement &me = _elements[getIndexFromIndices(ix, iy, 0)];
			LOWMSG2(int(me._subDomainLevel));
		}
		LOWMSG2("  " << _xlines[ix] << ":" << _xlines[ix+1] << "\n");
	}
	LOWMSG("");
	for(unsigned iz=0; iz<_zlines.size()-1; ++iz)
	{
		for(unsigned iy=0; iy<_ylines.size()-1;++iy)
		{
			const MeshElement &me = _elements[getIndexFromIndices(0, iy, iz)];
			LOWMSG2( char('A' + me._subDomainIdx));
		}
		LOWMSG2(" - ");
		for(unsigned iy=0; iy<_ylines.size()-1;++iy)
		{
			const MeshElement &me = _elements[getIndexFromIndices(0, iy, iz)];
			LOWMSG2(int(me._subDomainLevel));
		}
		LOWMSG2("  " << _zlines[iz] << ":" << _zlines[iz+1] << "\n");
	}
}

void Mesh::check() const
{

}

void Mesh::insert(LatticeSite *pLS, MeshElement *pEle)
{
	unsigned idx = (pEle? pEle->getIndex() : getIndexFromCoordinates(pLS->getCoordinates()));
	MeshElement &me = _elements[idx];
	pLS->_pElement = &me;
	pLS->_prev = 0;
	if(me._firstLS == 0)
	{
		me._firstLS = pLS;
		pLS->_next = 0;
	}
	else
	{
		me._firstLS->_prev = pLS;
		pLS->_next = me._firstLS;
		me._firstLS = pLS;
	}
}

void Mesh::remove(LatticeSite *pLS)
{
    MeshElement *pME = pLS->getElement();
    assert(pME->_firstLS);
    if(pLS->_prev == 0)
    {
        pME->_firstLS = pLS->_next;
        if(pLS->_next)
        	pLS->_next->_prev = 0;
    }
    else
    {
    	if(pLS->_next)
    		pLS->_next->_prev = pLS->_prev;
        pLS->_prev->_next = pLS->_next;
    }
}

void Mesh::insert(OKMC::Particle *pPart, MeshElement *pEle)
{
	unsigned idx = (pEle? pEle->getIndex() : getIndexFromCoordinates(pPart->getCoordinates()));
	if (idx >= _elements.size())
		ERRORMSG("Mesh::insert: Index overflow");

	trackLHF(idx, +1);
	MeshElement &me = _elements[idx];
	pPart->_pElement = &me;
	me._howManyParts++;
	pPart->_prev = 0;
	if(me._firstPart == 0)
	{
		me._firstPart = pPart;
		pPart->_next = 0;
	}
	else
	{
		me._firstPart->_prev = pPart;
		pPart->_next = me._firstPart;
		me._firstPart = pPart;
	}
	assert(pPart->_next != pPart);
}

void Mesh::remove(OKMC::Particle *pPart)
{
    MeshElement *pME = pPart->getElement();
    pME->_howManyParts--;
    trackLHF(pME->getIndex(), -1);
    assert(pME->_firstPart);
    if(pPart->_prev == 0)
    {
        pME->_firstPart = pPart->_next;
        if(pPart->_next)
        	pPart->_next->_prev = 0;
    }
    else
    {
    	if(pPart->_next)
    		pPart->_next->_prev = pPart->_prev;
    	pPart->_prev->_next = pPart->_next;
    }
}

unsigned Mesh::longHopFactor(unsigned idx) const
{
	if(_pDomain->_pMePar->_bLHF)
	{
		unsigned i,j,k;
		getIndicesFromIndex(idx, i, j, k);
		if(_xParticles[i] == 1 || _yParticles[j] == 1)
			return _longHopFactor[_elements[idx].getMaterial()].at(idx / _zCells);
	}
	return 1;
}

void Mesh::trackLHF(unsigned idx, int pm)
{
	unsigned ix,iy,iz;
	getIndicesFromIndex(idx, ix, iy, iz);

	_xParticles[ix]+=pm;
	if(ix != 0)
		_xParticles[ix-1]+=pm;
	else if(_periodicX)
		_xParticles[_xlines.size()-2]+=pm;
	if(ix < _xlines.size()-2)
		_xParticles[ix+1]+=pm;
	else if(_periodicX)
		_xParticles[0]+=pm;

	_yParticles[iy]+=pm;
	if(iy != 0)
		_yParticles[iy-1]+=pm;
	else if(_periodicY)
		_yParticles[_ylines.size()-2]+=pm;
	if(iy < _ylines.size()-2)
		_yParticles[iy+1]+=pm;
	else if(_periodicY)
		_yParticles[0]+=pm;
}

unsigned Mesh::getxIndex(float cx) const
{
	unsigned xInf = 0, xSup = _xlines.size() -1;
	unsigned nx;
	do
	{
	    nx = (xInf + xSup) / 2;
	    if(cx >= _xlines[nx+1])
	        xInf = nx;
	    else if(cx < _xlines[nx])
	        xSup = nx;
	    else
	        break;
	} while(true);
	return nx;
}

unsigned Mesh::getyIndex(float cy) const
{
	unsigned yInf = 0, ySup = _ylines.size() -1;
	unsigned ny;
	do
	{
	    ny = (yInf + ySup) / 2;
	    if(cy >= _ylines[ny+1])
	        yInf = ny;
	    else if(cy < _ylines[ny])
	        ySup = ny;
	    else
	        break;
	} while(true);
	return ny;
}

unsigned Mesh::getzIndex(float cz) const
{
	unsigned zInf = 0, zSup = _zlines.size() -1;
	unsigned nz;
	do
	{
	    nz = (zInf + zSup) / 2;
	    if(cz >= _zlines[nz+1])
	        zInf = nz;
	    else if(cz < _zlines[nz])
	        zSup = nz;
	    else
	        break;
	} while(true);
	return nz;
}

void Mesh::getIndicesFromIndex(unsigned idx, unsigned &ix, unsigned &iy, unsigned &iz) const
{
	ix = idx / _yzCells;
	unsigned rest = idx % _yzCells;
	iy = rest / _zCells;
	iz = rest % _zCells;
}

unsigned Mesh::getIndexFromCoordinates(const Coordinates &c) const
{
	int nx = getxIndex(c._x);
	int ny = getyIndex(c._y);
	int nz = getzIndex(c._z);
	return getIndexFromIndices(nx, ny, nz);
}

//takes BC into account
bool Mesh::fillLatticeNeighbors(const Coordinates &c, float dist, LSNeiList &neiList)
{
	const float dist2 = dist*dist;
	vector<MeshElement *> elems;
	fillNeighborsAllMat(c, dist, elems);
	for(vector<MeshElement *>::iterator it=elems.begin(); it!=elems.end(); ++it)
	{
		LatticeSite *pLS = (*it)->_firstLS;
		while(pLS)
		{
			Coordinates n = pLS->getCoordinates();
			float xdist = std::fabs(c._x - n._x);
			if(_periodicX && xdist > _xsize/2.) xdist -= _xsize;
			float ydist = std::fabs(c._y - n._y);
			if(_periodicY && ydist > _ysize/2.) ydist -= _ysize;
			float zdist = std::fabs(c._z - n._z);
			if(_periodicZ && zdist > _zsize/2.) zdist -= _zsize;
			float cdist2 = xdist*xdist + ydist*ydist + zdist*zdist;
			if(cdist2 < dist2)
				neiList.push_back(LSNeiInfo(pLS, cdist2));
			pLS = pLS->_next;
		}
	}
	return neiList.size();
}

void Mesh::getCenter(Coordinates &c, unsigned idx) const
{
	unsigned ix, iy, iz;
	getIndicesFromIndex(idx, ix, iy, iz);
	c._x = (_xlines[ix] + _xlines[ix+1])/2.;
	c._y = (_ylines[iy] + _ylines[iy+1])/2.;
	c._z = (_zlines[iz] + _zlines[iz+1])/2.;
}

void Mesh::getCorners(unsigned idx, Coordinates &m, Coordinates &M) const
{
	unsigned ix, iy, iz;
	getIndicesFromIndex(idx, ix, iy, iz);
	m._x = _xlines[ix]; M._x = _xlines[ix+1];
	m._y = _ylines[iy]; M._y = _ylines[iy+1];
	m._z = _zlines[iz]; M._z = _zlines[iz+1];
}

double Mesh::getVolume(unsigned idx) const
{
	unsigned ix, iy, iz;
	getIndicesFromIndex(idx, ix, iy, iz);
	return (_xlines[ix+1]  - _xlines[ix]) * (_ylines[iy+1] - _ylines[iy]) * (_zlines[iz+1] - _zlines[iz]);
}

void Mesh::getRandomPosition(unsigned idx, Coordinates &m, float r1, float r2, float r3) const
{
	unsigned ix, iy, iz;
	getIndicesFromIndex(idx, ix, iy, iz);
	m._x = _xlines[ix] + r1*(_xlines[ix+1] - _xlines[ix]);
	m._y = _ylines[iy] + r2*(_ylines[iy+1] - _ylines[iy]);
	m._z = _zlines[iz] + r3*(_zlines[iz+1] - _zlines[iz]);
}

//returns REJECTED or OK
//Applies Boundary Conditions (BC)
//changes "origin" when applying BC
Mesh::JUMP_ACTIONS Mesh::jumpPosition(Coordinates &c, const Coordinates &lambda, MeshElement *&pEle,
		const std::pair<P_TYPE, unsigned> *jumpType, Coordinates &o, unsigned longHopF)
{
	Coordinates old = c;
	unsigned ix, iy, iz;
	getIndicesFromIndex(pEle->getIndex(), ix, iy, iz);
	if(jumpType && jumpType->first != Kernel::UNDEFINED_TYPE)
	{
		std::map<std::pair<P_TYPE, unsigned>, unsigned>::iterator it = pEle->_hops.find(*jumpType);
		if(it == pEle->_hops.end())
			pEle->_hops[*jumpType] = longHopF;
		else
			it->second+=longHopF;
	}

	c += lambda;
	bool bElementChanged = false;
	if(c._x < _xlines[ix])
	{
		if(c._x < _xlines[0])
		{
			if(_periodicX)
			{
				ix = _xlines.size() -2;
				c._x += _xsize;
				o._x += _xsize;
				if(c._x == _xlines.back())
					c._x -= 1e-5; //FP issues
			}
			else
				goto rejected;
		}
		while(c._x < _xlines[ix]) //update ix.
			ix--;
		bElementChanged = true;
	}
	else if(c._x >= _xlines[ix+1])
	{
		if(c._x >= _xlines.back())
		{
			if(_periodicX)
			{
				ix = 0;
				c._x -= _xsize;
				o._x -= _xsize;
			}
			else
				goto rejected;
		}
		while(c._x >= _xlines[ix+1])
			ix++;
		bElementChanged = true;
	}

	if(c._y < _ylines[iy])
	{
		if(c._y < _ylines[0])
		{
			if(_periodicY)
			{
				iy = _ylines.size() -2;
				c._y += _ysize;
				o._y += _ysize;
				if(c._y == _ylines.back())
					c._y -= 1e-5; //FP issues
			}
			else
				goto rejected;
		}
		while(c._y < _ylines[iy])
			iy--;
		bElementChanged = true;
	}
	else if(c._y >= _ylines[iy+1])
	{
		if(c._y >= _ylines.back())
		{
			if(_periodicY)
			{
				iy = 0;
				c._y -= _ysize;
				o._y -= _ysize;
			}
			else
				goto rejected;
		}
		while(c._y >= _ylines[iy+1])
			iy++;
		bElementChanged = true;
	}

	if(c._z < _zlines[iz])
	{
		if(c._z < _zlines[0])
		{
			if(_periodicZ)
			{
				iz = _zlines.size() - 2;
				c._z += _zsize;
				o._z += _zsize;
				if(c._z == _zlines.back())
					c._z -= 1e-5; //FP issues
			}
			else
				goto rejected;
		}
		while(c._z < _zlines[iz])
			iz--;
		bElementChanged = true;
	}
	else if(c._z >= _zlines[iz+1])
	{
		if(c._z >= _zlines.back())
		{
			if(_periodicZ)
			{
				iz = 0;
				c._z -= _zsize;
				o._z -= _zsize;
			}
			else
				goto rejected;
		}
		while(c._z >= _zlines[iz+1])
			iz++;
		bElementChanged = true;
	}
	if(bElementChanged)
		pEle = getElement(getIndexFromIndices(ix, iy, iz));
	return JUMP_OK;

rejected:
	c = old;
	return JUMP_REJECTED;
}

//Takes BC into account
void Mesh::middle(Coordinates &c, const Coordinates &c1, const Coordinates &c2) const
{
	c = c2;
	c+=c1;
	c*=.5;
	if(_periodicX && std::fabs(c2._x - c1._x) >= _xsize/2.)
	{
		c._x += _xsize/2.;
		if(c._x >= _xlines.back())
			c._x -= _xsize;
	}
	if(_periodicY && std::fabs(c2._y - c1._y) >= _ysize/2.) //periodic in Y
	{
		c._y += _ysize/2.;
		if(c._y >= _ylines.back())
			c._y -= _ysize;
	}
	if(_periodicZ && std::fabs(c2._z - c1._z) >= _zsize/2.) //periodic in Z
	{
		c._z += _zsize/2.;
		if(c._z >= _zlines.back())
			c._z -= _zsize;
	}
}

//it takes periodic bc into account
void Mesh::fillNeighborsOneMat(const MeshElement *pEle, const Coordinates &c, float radius, vector<MeshElement *> &elems)
{
	M_TYPE mt = pEle->getMaterial();
	unsigned ix, iy, iz;
	getIndicesFromIndex(pEle->getIndex(), ix, iy, iz);

	//obtain the neighboring elements
	unsigned maxX = ix+1, minX = ix;
	float limit = c._x - radius;
	bool overwrite = false;
	while(limit < _xlines[minX])
	{
		minX--;
		if(int(minX) < 0)
		{
			if(_periodicX)
			{
				minX  = _xlines.size() -1;
				limit += _xsize;
			}
			else
			{
				minX = 0;
				break;
			}
		}
		if(minX == maxX)
			overwrite = true;
	}
	limit = c._x + radius;
	while(limit >= _xlines[maxX])
	{
		maxX++;
		if(maxX >= _xlines.size())
		{
			if(_periodicX)
			{
				maxX = 1;
				limit -= _xsize;
			}
			else
			{
				maxX = _xlines.size()-1;
				break;
			}
		}
		if(minX == maxX)
			overwrite = true;
	}
	if(overwrite)
	{
		WARNINGMSG("Defect large enough to see itself using PBC in x");
		minX = 0;
		maxX = _xlines.size() - 1;
	}
	// Y
	overwrite = false;
	unsigned maxY = iy+1, minY = iy;
	limit = c._y - radius;
	while(limit < _ylines[minY])
	{
		minY--;
		if(int(minY) < 0)
		{
			if(_periodicY)
			{
				minY  = _ylines.size() -1;
				limit += _ysize;
			}
			else
			{
				minY = 0;
				break;
			}
		}
		if(minY == maxY)
			overwrite = true;
	}
	limit = c._y + radius;
	while(limit >= _ylines[maxY])
	{
		maxY++;
		if(maxY >= _ylines.size())
		{
			if(_periodicY)
			{
				maxY = 1;
				limit -= _ysize;
			}
			else
			{
				maxY = _ylines.size()-1;
				break;
			}
		}
		if(minY == maxY)
			overwrite = true;
	}
	if(overwrite)
	{
		WARNINGMSG("Defect large enough to see itself using PBC in y");
		minY = 0;
		maxY = _ylines.size() - 1;
	}
	// Z
	overwrite = false;
	unsigned maxZ = iz+1, minZ = iz;
	limit = c._z - radius;
	while(limit < _zlines[minZ])
	{
		minZ--;
		if(int(minZ) < 0)
		{
			if(_periodicZ)
			{
				minZ  = _zlines.size() -1;
				limit += _zsize;
			}
			else
			{
				minZ = 0;
				break;
			}
		}
		if(minZ == maxZ)
			overwrite = true;
	}
	limit = c._z + radius;
	while(limit >= _zlines[maxZ])
	{
		maxZ++;
		if(maxZ >= _zlines.size())
		{
			if(_periodicZ)
			{
				maxZ = 1;
				limit -= _zsize;
			}
			else
			{
				maxZ = _zlines.size()-1;
				break;
			}
		}
		if(minZ == maxZ)
			overwrite = true;
	}
	if(overwrite)
	{
		WARNINGMSG("Defect large enough to see itself using PBC in z");
		minZ = 0;
		maxZ = _zlines.size() - 1;
	}
	elems.reserve(27);
	for(unsigned indx = minX; indx != maxX; ++indx)
	{
		if(indx == _xlines.size()-1)
			indx=0;
		for(unsigned indy = minY; indy != maxY; ++indy)
		{
			if(indy == _ylines.size()-1)
				indy=0;
			for(unsigned indz = minZ; indz != maxZ; ++indz)
			{
				if(indz == _zlines.size()-1)
					indz=0;
				unsigned index = getIndexFromIndices(indx, indy, indz);
				if(_elements[index].getMaterial() == mt)
					elems.push_back(&_elements[index]);
			}
		}
	}
}

//fills neighbors for all materials around c, distance dist
void Mesh::fillNeighborsAllMat(const Coordinates &c, float dist, vector<MeshElement *> &elems)
{
	float newdist = c[0]-dist;
	if(newdist < _xlines.front())
	{
		newdist = (_periodicX? _xsize + newdist : _xlines.front());
		if(newdist >= _xlines.back()) newdist = _xlines.front();
	}
	const unsigned ix = getxIndex(newdist);
	newdist = c[0]+dist;
	if(newdist >= _xlines.back())
	{
		newdist = (_periodicX? newdist - _xsize : _xlines.back() - 1e-4);
		if(newdist < _xlines.front()) newdist = _xlines.back() - 1e-4;
	}
	const unsigned Ix = getxIndex(newdist);
	newdist = c[1]-dist;
	if(newdist < _ylines.front())
	{
		newdist = (_periodicY? _ysize + newdist : _ylines.front());
		if(newdist >= _ylines.back()) newdist = _ylines.front(); //to avoid floating issues when applying BC
	}
	const unsigned iy = getyIndex(newdist);
	newdist = c[1]+dist;
	if(newdist >= _ylines.back())
	{
		newdist = (_periodicY? newdist - _ysize : _ylines.back() - 1e-4);
		if(newdist < _ylines.front()) newdist = _ylines.back() - 1e-4;
	}
	const unsigned Iy = getyIndex(newdist);
	newdist = c[2]-dist;
	if(newdist < _zlines.front())
	{
		newdist = (_periodicZ? _zsize + newdist : _zlines.front());
		if(newdist >= _zlines.back()) newdist = _zlines.front();
	}
	const unsigned iz = getzIndex(newdist);
	newdist = c[2]+dist;
	if(newdist >= _zlines.back())
	{
		newdist = (_periodicZ? newdist - _zsize : _zlines.back() - 1e-4);
		if(newdist < _zlines.front()) newdist = _zlines.back() - 1e-4;
	}
	const unsigned Iz = getzIndex(newdist);

	for(unsigned i=ix; (ix <= Ix? (i<=Ix) : (i >= ix || i <= Ix) ); ++i)
	{
		if(i == _xlines.size() -1) i=0;
		for(unsigned j=iy; (iy <= Iy? (j <=Iy) : (j >= iy || j <= Iy) ); ++j)
		{
			if(j == _ylines.size() -1) j = 0;
			for(unsigned k=iz; (iz <= Iz? (k <= Iz) : (k >=iz || k <= Iz) ); ++k)
			{
				if(k == _zlines.size() - 1) k = 0;
				elems.push_back(&_elements[getIndexFromIndices(i, j, k)]);
			}
		}
	}
}

void Mesh::fillInterfaces(const Coordinates &c, float dist, vector<OKMC::Interface *> &ints)
{
	vector<MeshElement *> elems;
	fillNeighborsAllMat(c, dist, elems);
	for(vector<MeshElement *>::iterator it = elems.begin(); it != elems.end(); ++it)
	{
		const vector<OKMC::Interface *> interfaces = (*it)->getInterfaces();
		for(vector<OKMC::Interface *>::const_iterator itInt = interfaces.begin(); itInt != interfaces.end(); ++itInt)
			if(std::find(ints.begin(), ints.end(), *itInt) == ints.end())
				ints.push_back(*itInt);
	}

}


void Mesh::fillAdjacentNeighborsOneMat(const MeshElement *pEle, vector<MeshElement *> &neighs)
{
	unsigned x, y, z;
	getIndicesFromIndex(pEle->getIndex(), x, y, z);
	for (signed i = -1; i <= 1 ; i++)
		for (signed j = -1; j <= 1; j++)
			for (signed k = -1; k <= 1; k++)
			{
				MeshElement * ele;
				if(step(x, y, z, i,  j, k) == NULL || (i == 0 && j == 0 && k == 0)) continue;
				ele = const_cast<MeshElement*>(step(x, y, z, i, j, k));                                
				if(ele->getMaterial() == pEle->getMaterial())
					neighs.push_back(ele);
			}
}
// Used for smoothing, returns the elements in the sides, edges and corners of the given one
// Only returns elements of the same material of the central element
// Takes into account BC
void Mesh::fillNeighborsTopology(MeshElement *pEle, vector<const MeshElement *> &sides,
		 vector<const MeshElement *> &edges, vector<const MeshElement *> &corners) const
{
	unsigned x, y, z;
	getIndicesFromIndex(pEle->getIndex(), x, y, z);
	for (signed i = -1; i <= 1 ; i++)
		for (signed j = -1; j <= 1; j++)
			for (signed k = -1; k <= 1; k++)
			{
				const MeshElement * ele;
				if(step(x, y, z, i,  j, k) == NULL || (i == 0 && j == 0 && k == 0)) continue;
				ele = step(x, y, z, i, j, k);
				if((i == 0 && j == 0) || (j == 0 && k == 0) || (i == 0 && k == 0))
				{
					if(ele->getMaterial() == pEle->getMaterial())
						sides.push_back(ele);
				}
				else if(i != 0 && j != 0 && k != 0)
				{
					if(ele->getMaterial() == pEle->getMaterial())
						corners.push_back(ele);
				}
				else
				{
					if(ele->getMaterial() == pEle->getMaterial())
						edges.push_back(ele);
				}
			}
}

// Returns elements of different material than the central element also
// Takes into account BC
// if mt or basicMt are different from MAX_MATERIALS, discard boxes without that material or basic material
void Mesh::fillAdjacentNeighbors(const MeshElement *pEle, vector<MeshElement *> &all, M_TYPE mt, M_TYPE basicMt)
{
	const IO::ParameterManager *pPM = Domains::global()->PM();
	unsigned x, y, z;
	getIndicesFromIndex(pEle->getIndex(), x, y, z);
	for (signed i = -1; i <= 1 ; i++)
		for (signed j = -1; j <= 1; j++)
			for (signed k = -1; k <= 1; k++)
			{
				MeshElement * ele;
				if(!step(x, y, z, i,  j, k) || (i == 0 && j == 0 && k == 0)) continue;
				ele = const_cast<MeshElement*>(step(x, y, z, i, j, k));
				if(mt != MAX_MATERIALS && ele->getMaterial() != mt) continue;
				if(basicMt != MAX_MATERIALS && pPM->getMaterial(ele->getMaterial())._basicMaterial != basicMt) continue;
				all.push_back(ele);
			}
}


//it takes BC into account
//dist is automatically computed
void Mesh::getInteractions(SubDomain *pSub, const OKMC::Particle *part, vector<OKMC::Particle *> &parts)
{
	vector<MeshElement *> elems;
	const Coordinates &c = part->getCoordinates();
	float radius = Domains::global()->PM()->getMaximumCaptureRadius(part->getElement()->getMaterial());
	fillNeighborsOneMat(part->getElement(), c, radius, elems);

	//compute distances and insert right ones.
	for(vector<MeshElement *>::iterator it=elems.begin(); it!=elems.end(); ++it)
	{
		OKMC::Particle *ppart = (*it)->_firstPart;
		while(ppart)
		{
			if(ppart->getDefect() != part->getDefect())
				if(part->getDefect()->canInteract(pSub, ppart, part))
					parts.push_back(ppart);
			ppart = ppart->_next;
		}
	}
}

//accepts BC
void Mesh::fillNeighborsOneMat(const MeshElement *pEle, const Coordinates &c, float dist, vector<OKMC::Particle *> &parts)
{
	vector<MeshElement *> elems;
	float radius = Domains::global()->PM()->getMaximumCaptureRadius(pEle->getMaterial());
	fillNeighborsOneMat(pEle, c, radius, elems);

	//compute distances and insert right ones.
	const float dist2 = dist*dist;
	for(vector<MeshElement *>::iterator it=elems.begin(); it!=elems.end(); ++it)
	{
		OKMC::Particle *ppart = (*it)->_firstPart;
		while(ppart)
		{
			Coordinates n = ppart->getCoordinates();
			float xdist = std::fabs(c._x - n._x);
			if(_periodicX && xdist > _xsize/2.) xdist -= _xsize;
			float ydist = std::fabs(c._y - n._y);
			if(_periodicY && ydist > _ysize/2.) ydist -= _ysize;
			float zdist = std::fabs(c._z - n._z);
			if(_periodicZ && zdist > _zsize/2.) zdist -= _zsize;
			if( xdist*xdist + ydist*ydist + zdist*zdist < dist2)
				parts.push_back(ppart);
			ppart = ppart->_next;
		}
	}
}

/*
 * Returns n-c with periodic BC
 * writes the solution in n (overwrites n)
 */
void Mesh::setPeriodicRelative(const Coordinates &c, Coordinates &n) const
{
	n._x -= c._x;
	if(_periodicX && n._x > _xsize/2.)
		n._x -= _xsize;
	else if(_periodicX && n._x < -_xsize/2.)
		n._x += _xsize;
	n._y -= c._y;
	if(_periodicY && n._y > _ysize/2.)
		n._y -= _ysize;
	else if(_periodicY && n._y < -_ysize/2.)
		n._y += _ysize;
	n._z -= c._z;
	if(_periodicZ && n._z > _zsize/2.)
		n._z -= _zsize;
	else if(_periodicZ && n._z < -_zsize/2.)
		n._z += _zsize;
}

/* Boost version */
void Mesh::setPeriodicRelative(const ublas::vector<float> &c, ublas::vector<float> &n) const
{
	n -= c;
	if(_periodicX && n(0) > _xsize/2.)
		n(0) -= _xsize;
	else if(_periodicX && n(0) < -_xsize/2.)
		n(0) += _xsize;
	if(_periodicY && n(1) > _ysize/2.)
		n(1) -= _ysize;
	else if(_periodicY && n(1) < -_ysize/2.)
		n(1) += _ysize;
	if(_periodicZ && n(2) > _zsize/2.)
		n(2) -= _zsize;
	else if(_periodicZ && n(2) < -_zsize/2.)
		n(2) += _zsize;
}

void Mesh::setPeriodicRelative(const ublas::vector<double> &c, ublas::vector<double> &n) const
{
	n -= c;
	if(_periodicX && n(0) > _xsize/2.)
		n(0) -= _xsize;
	else if(_periodicX && n(0) < -_xsize/2.)
		n(0) += _xsize;
	if(_periodicY && n(1) > _ysize/2.)
		n(1) -= _ysize;
	else if(_periodicY && n(1) < -_ysize/2.)
		n(1) += _ysize;
	if(_periodicZ && n(2) > _zsize/2.)
		n(2) -= _zsize;
	else if(_periodicZ && n(2) < -_zsize/2.)
		n(2) += _zsize;
}

void Mesh::getDomain(Coordinates &m, Coordinates &M) const
{
	m._x = _xlines.front(); M._x = _xlines.back();
	m._y = _ylines.front(); M._y = _ylines.back();
	m._z = _zlines.front(); M._z = _zlines.back();
}

bool Mesh::isInDomain(const Coordinates &c) const
{
	if(c._x < _xlines.front() || c._x >= _xlines.back() ||
	   c._y < _ylines.front() || c._y >= _ylines.back() ||
	   c._z < _zlines.front() || c._z >= _zlines.back())
			return false;
	return true;
}

LatticeSite * Mesh::findLS(const Coordinates &c, float dist)
{
	LSNeiList neiList;
	fillLatticeNeighbors(c, dist, neiList);
	assert(neiList.size() <= 1);
	if(neiList.size())
		return neiList.back()._pLS;
	else
		return 0;
}

//applies per. boun. cond. and sees if there are potential interfaces. Updates pEle
//requires a pEle pointing the to right material, but nothing else.
Mesh::JUMP_ACTIONS Mesh::checkMove(MeshElement *&pEle, Coordinates &to)
{
	M_TYPE mt1 = pEle->getMaterial();
	while(to._x < _xlines.front() || to._x >= _xlines.back())
	{
		if(!_periodicX)
			return JUMP_REJECTED;
		to._x += (to._x <= _xlines.front()? _xsize : -_xsize);
	}
	while(to._y < _ylines.front() || to._y >= _ylines.back())
	{
		if(!_periodicY)
			return JUMP_REJECTED;
		to._y += (to._y <= _ylines.front()? _ysize : -_ysize);
	}
	while(to._z < _zlines.front() || to._z >= _zlines.back())
	{
		if(!_periodicZ)
			return JUMP_REJECTED;
		to._z += (to._z <= _zlines.front()? _zsize : -_zsize);
	}
	pEle = getElement(getIndexFromCoordinates(to));
	if(pEle->getMaterial() == mt1)
		return JUMP_OK;
	return JUMP_INTERFACE;
}

//changes a material and updates the interfaces
void Mesh::changeMaterial(Kernel::SubDomain *pSub, unsigned idx, M_TYPE to)
{
	MeshElement *pME = getElement(idx);
	M_TYPE from = pME->_mat;

	if(pME->getInterfaces().size())
	{
		vector<OKMC::Particle *> part;
		unsigned ix, iy, iz;
		getIndicesFromIndex(idx, ix, iy, iz);
		unsigned nx  = _pDomain->_pMesh->getLines(0).size() - 2;
		unsigned ny  = _pDomain->_pMesh->getLines(1).size() - 2;
		unsigned nz  = _pDomain->_pMesh->getLines(2).size() - 2;

		while (pME->_interfaces.size())
		{
			vector<OKMC::Particle *> tmp = pME->_interfaces[0]->getParticles();
			for (unsigned j = 0; j < tmp.size(); ++j)
				// a copy of the particle is necessary because it will be deleted in interface destructor
				part.push_back(new OKMC::Particle(tmp[j]->getPType(), tmp[j]->getCoordinates(), tmp[j]->getDefect(), tmp[j]->getOrig()));
			_pDomain->_pRM->remove(pME->_interfaces[0], pME->_interfaces[0]->getElement(0));
			delete pME->_interfaces[0];
		}
		pME->_mat = to; //we change the material here, because we are creating interfaces. We'll change it again just after we finish
		if (ix >= 1)
		{
			Kernel::MeshElement *pME2 = getElement(getIndexFromIndices(ix - 1, iy, iz));                        
			if (pME2->getMaterial() != to)
			{
				OKMC::Interface *pIF = new OKMC::Interface(_pDomain, pME2, pME, 0);
				_pDomain->_pRM->insert(pIF, pME2);
			}
		}
		if (iy >= 1)
		{
			Kernel::MeshElement *pME2 = getElement(getIndexFromIndices(ix, iy - 1, iz));
			if (pME2->getMaterial() != to)
			{
				OKMC::Interface *pIF = new OKMC::Interface(_pDomain, pME2, pME, 1);
				_pDomain->_pRM->insert(pIF, pME2);
			}
		}
		if (iz >= 1)
		{
			Kernel::MeshElement *pME2 = getElement(getIndexFromIndices(ix, iy, iz - 1));
			if (pME2->getMaterial() != to)
			{
				OKMC::Interface *pIF = new OKMC::Interface(_pDomain, pME2, pME, 2);
				_pDomain->_pRM->insert(pIF, pME2);
			}
		}
		if (ix < nx)
		{
			Kernel::MeshElement *pME2 = getElement(getIndexFromIndices(ix + 1, iy, iz));
			if (pME2->getMaterial() != to)
			{
				OKMC::Interface *pIF = new OKMC::Interface(_pDomain, pME, pME2, 0);
				_pDomain->_pRM->insert(pIF, pME);
			}
		}
		if (iy < ny)
		{
			Kernel::MeshElement *pME2 = getElement(getIndexFromIndices(ix, iy + 1, iz));
			if (pME2->getMaterial() != to)
			{
				OKMC::Interface *pIF = new OKMC::Interface(_pDomain, pME, pME2, 1);
				_pDomain->_pRM->insert(pIF, pME);
			}
		}
		if (iz < nz)
		{
			Kernel::MeshElement *pME2 = getElement(getIndexFromIndices(ix, iy, iz + 1));
			if (pME2->getMaterial() != to)
			{
				OKMC::Interface *pIF = new OKMC::Interface(_pDomain, pME, pME2, 2);
				_pDomain->_pRM->insert(pIF, pME);
			}
		}

		if (pME->_interfaces.size())
		{
			for (unsigned i = 0; i < part.size(); ++i)
			{
				int r = (pSub->_rng.rand() * pME->_interfaces.size());
				pME->_interfaces[r]->insertParticle(part[i]);
				// we should check that the particle is inserted at an interface of the same type?!
				// in that case we have to take care because an interface of the same type might not exist
			}

			for (unsigned i = 0; i < pME->_interfaces.size(); ++i)
				_pDomain->_pRM->update(pME->_interfaces[i], pME->_interfaces[i]->getElement(0));

		}
	}
	pME->_mat = from;
	//map particles from one material to the other
	const IO::MaterialProperties &fromProp = Domains::global()->PM()->getMaterial(from);
	const IO::MaterialProperties &  toProp = Domains::global()->PM()->getMaterial(to);

	vector<OKMC::Particle *> partSet;  //cannot be a set because that will introduce
		//a dependency with the pointer number, and thus the same binary will produce
		//different results even with the same seed.
	for (OKMC::Particle *p = pME->_firstPart; p; p = p->getNext())
		partSet.push_back(p);

	// solid to gas (vaporization)            ||  //crystalline to amorphous (amorphization)
	pME->_mat = to;

	if(Domains::global()->PM()->isGasLike(to) || (fromProp._bAmorphous == false && toProp._bAmorphous == true)
			|| (Domains::global()->PM()->isGasLike(from) && toProp._bEpitaxy) ) //gas to other (epitaxy)
	{
		transfer(pSub, partSet, pME, from, to);
		return;
	}

	// amorphous to crystalline (recrystallization)
	if(fromProp._bAmorphous == true && toProp._bAmorphous == false)
	{
		//the model assumes that the interface "brooms" all particles
		for(vector<OKMC::Particle *>::iterator it=partSet.begin(); it != partSet.end(); ++it)
		{
			int r = pSub->_rng.rand() * pME->_interfaces.size();
			if(pME->_interfaces.size() && pME->_interfaces[r]->canInteract(pSub, *it, 0))
				pME->_interfaces[r]->interact(pSub, (*it)->getDefect(), 0, partSet);
		}
		//rest of particles
		transfer(pSub, partSet, pME, from, to);
		return;
	}

	transfer(pSub, partSet, pME, from, to);
	MEDMSG(fromProp._name << "->" << toProp._name);
}

//transfers the particle between materials when the material changes.
void Mesh::transfer(SubDomain *pSub, vector<OKMC::Particle *> &partSet, MeshElement *pME, M_TYPE from, M_TYPE to)
{
	while(partSet.size())
	{
		P_TYPE transferType = partSet.back()->getPType();
		if( Domains::global()->PM()->getPPos(transferType) == NO_POS ||  // Mobile particle cannot be NOPOS
			!Domains::global()->PM()->isParticleDefined(transferType, to))
		{
			P_TYPE fam = Domains::global()->PM()->getFamily(transferType);
			transferType = Domains::global()->PM()->getParticle(to, fam, POS_0);
		}
		if(!Domains::global()->PM()->isParticleDefined(transferType, to))
			transferType = UNDEFINED_TYPE;

		Kernel::Coordinates c = partSet.back()->getCoordinates(), orig = partSet.back()->getOrig();
		partSet.back()->getDefect()->deletePart(pSub, partSet); //takes care of partSet
		if(transferType != UNDEFINED_TYPE)
			new OKMC::MobileParticle(pSub, transferType, 0, pME->getDomain(), c, orig);
	}
}

void Mesh::resetDiffusivity()
{
	//reset orig positions for diffusivity calculation.
	for(iterator it=begin(); it!=end(); ++it)
	{
		OKMC::Particle *part = it->getFirstPart();
		while(part)
		{
			part->setOrig(part->getCoordinates());
			part = part->getNext();
		}
	}
}

const MeshElement * Mesh::step(unsigned ix, unsigned iy, unsigned iz, int xStep, int yStep, int zStep) const
{
	if(xStep == 1)
	{
		if(_periodicX)
			ix < _xlines.size() - 2 ? ix++ : ix = 0;
		else
		{
			if(ix < _xlines.size() - 2)
				ix++;
			else
				return NULL;
		}
	}
	if(yStep == 1)
	{
		if(_periodicY)
			iy < _ylines.size() - 2 ? iy++ : iy = 0;
		else
		{
			if(iy < _ylines.size() - 2)
				iy++;
			else
				return NULL;
		}
	}
	if(zStep == 1)
	{
		if(_periodicZ)
			iz < _zlines.size() - 2 ? iz++ : iz = 0;
		else
		{
			if(iz < _zlines.size() - 2)
				iz++;
			else
				return NULL;
		}
	}
	if(xStep == -1)
	{
		if(_periodicX)
			ix > 0 ? ix-- : ix = _xlines.size() - 2;
		else
		{
			if(ix > 0)
				ix--;
			else
				return NULL;
		}
	}
	if(yStep == -1)
	{
		if(_periodicY)
			iy > 0 ? iy-- : iy = _ylines.size() - 2;
		else
		{
			if(iy > 0)
				iy--;
			else
				return NULL;
		}
	}
	if(zStep == -1)
	{
		if(_periodicZ)
			iz > 0 ? iz-- : iz = _zlines.size() - 2;
		else
		{
			if(iz > 0)
				iz--;
			else
				return NULL;
		}
	}
	return &_elements[getIndexFromIndices(ix, iy, iz)];
}

void Mesh::buildNeighbors()
{
	for(vector<MeshElement>::iterator it=_elements.begin(); it!=_elements.end(); ++it)
		fillNeighborsTopology(&(*it), it->_sides, it->_edges, it->_corners);
}

/* This function is called when a cluster containing impurities either produces a "bump" in the surface, or
 * touches such surface, or produces a "trench"
 *
 * If a "bump" is produced, then howMany is +1
 * If a hole is produced, then howMany is a negative number.
 *
 * The goal is to look for the closer free surface and operate with it.
 */

void Mesh::changeNumberOfCells(Kernel::SubDomain *pSub, MeshElement *pEle, int howMany)
{
	//checks for all gas cells belonging to an interface, and computes the closest one to the pEle
	double distance = 1e30; //absurd initial number
	vector<MeshElement>::iterator toOperate = _elements.end();
	M_TYPE gasMt = MAX_MATERIALS;
	for(vector<MeshElement>::iterator it=_elements.begin(); it!=_elements.end(); ++it)
	{
		bool cond = Domains::global()->PM()->isGasLike(it->getMaterial());
		if(cond)
			gasMt = it->getMaterial();
		if(howMany == -1) //needs a material element to contract
			cond = it->getMaterial() == pEle->getMaterial();
		if(howMany != -1 && howMany != +1)
			ERRORMSG("changeNumberOfCells called with improper argument!");
		if(cond && it->getInterfaces().size())
		{
			Coordinates c = pEle->getCoords(.5, .5, .5);
			setPeriodicRelative(it->getCoords(.5, .5, .5), c);
			double dist = c*c;
			if(c*c < distance)
			{
				distance = dist;
				toOperate = it;
			}
		}
	}
	assert(gasMt != MAX_MATERIALS);
	if(toOperate == _elements.end())  //something failed
		return;

	M_TYPE newMt = (howMany == +1? pEle->getMaterial() : gasMt);
	changeMaterial(pSub, toOperate->getIndex(), newMt);
}

// Moves a lattice atom B from the AB alloy from MeshElement from to MeshElement to with the given probability
// if not moves an A atom
void Mesh::jumpBorA(MeshElement * from, MeshElement * to, double p)
{
    if (_pDomain->_rng_dom.rand() < p) // B jumps from ME*from to ME*to
    {
        from->decBAtoms();
        to->incBAtoms();
    }
    else // A jumps from ME*from to ME*to
    {
        from->decAAtoms();
        to->incAAtoms();
    }
}
