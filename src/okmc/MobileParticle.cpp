/*
 * MobileParticle.cpp
 *
 *  Created on: Feb. 24, 2011
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
 *
 *
 *  Keeps diffusion for impurities only
 *
 */

#include "MobileParticle.h"
#include "kernel/Constants.h"
#include "kernel/Material.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/MeshElement.h"
#include "kernel/Mesh.h"
#include "kernel/RateManager.h"
#include "electrostatics/ParticleToNodeHandler.h"
#include "kernel/MeshParam.h"
#include "lkmc/Lattice.h"
#include "MobileParticleParam.h"
#include "AlloyParam.h"
#include "Interface.h"
#include "Cluster.h"
#include "io/ParameterManager.h"
#include "kernel/StateManager.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <stdlib.h>
#include "io/FileParameters.h"

using Kernel::Mesh;
using std::vector;
using Kernel::Coordinates;
using namespace Electrostatics;
using std::map;

using Kernel::Coordinates;
using Kernel::P_TYPE;
using Kernel::P_POS;
using Kernel::POS_0;
using Kernel::POS_1;
using Kernel::ID;

namespace OKMC {

MobileParticle::MobileParticle(Kernel::SubDomain *pSub, P_TYPE type,unsigned state, Kernel::Domain *p, const Kernel::Coordinates &c, const Kernel::Coordinates &o)
: Defect(p), Particle(type, c, this, o)
{
	assert(_ptype != Kernel::UNDEFINED_TYPE);
	_state=state;
	_pDomain->_pMesh->insert(this, 0);
	if (Domains::global()->IsPoisson())
		_pDomain->_pMesh->getParticleToNodeHandler()->insert(this);
	assert(Domains::global()->PM()->isParticleDefined(type, _pElement->getMaterial()));
	_pDomain->_pRM->insert(this, _pElement);
	setAxes(pSub);
	if (Domains::global()->PM()->getMaterial(_pElement->getMaterial())._selfDiffusion && 
            Domains::global()->PM()->getIorV(_pElement->getMaterial(), type) == Kernel::IV_I)
	{
		if(pSub->_rng.rand() > _pElement->getAlloyFraction())
		{
			_pElement->incAAtoms();
			pSub->_evLog.performed(_pElement->getMaterial(), Event::CELL, _ptype, 0, 0, _state);
			_pElement->_ABalance++;
		}
		else
		{
			_pElement->incBAtoms();
                        pSub->_evLog.performed(_pElement->getMaterial(), Event::CELL, _ptype, 0, 2, _state);
			_pElement->_BBalance++;
		}
	}
	if (Domains::global()->PM()->getMaterial(_pElement->getMaterial())._selfDiffusion && 
            Domains::global()->PM()->getIorV(_pElement->getMaterial(), type) == Kernel::IV_V)
	{
		if(pSub->_rng.rand() > _pElement->getAlloyFraction())
		{
			_pElement->decAAtoms();
                        pSub->_evLog.performed(_pElement->getMaterial(), Event::CELL, _ptype, 0, 1, _state);
			_pElement->_ABalance--;
		}
		else
		{
			_pElement->decBAtoms();
                        pSub->_evLog.performed(_pElement->getMaterial(), Event::CELL, _ptype, 0, 3, _state);
			_pElement->_BBalance--;
		}
	}
#ifdef NUMODEL
	for (unsigned i = 0; i < 6; i++)
		_AorB[i] = _pDomain->_rng_dom.rand() > _pElement->getAlloyFraction() ? true : false;
#endif
}

MobileParticle::MobileParticle(std::istream &is) : Defect(is), Particle(is, this)
{
	is >> _state;
	_pDomain->_pMesh->insert(this, 0);
	if (Domains::global()->IsPoisson())
		_pDomain->_pMesh->getParticleToNodeHandler()->insert(this);
	_pDomain->_pRM->insert(this, _pElement);
	setAxes(_pDomain->_pRM->getSubDomain(0)); //anyone should work
}

void MobileParticle::restart(std::ostream &os) const
{
	Defect::restart(os);
	Particle::restart(os);
	os << _state;
}

MobileParticle::~MobileParticle()
{        
        Kernel::SubDomain *pSub = _pDomain->_pRM->getSubDomain(_pElement->getSubDomainIdx());
    	if (Domains::global()->PM()->getMaterial(_pElement->getMaterial())._selfDiffusion && 
            Domains::global()->PM()->getIorV(_pElement->getMaterial(), _ptype) == Kernel::IV_I)
	{
		if(_pDomain->_rng_dom.rand() > _pElement->getAlloyFraction())
		{
			_pElement->decAAtoms();                        
			pSub->_evLog.performed(_pElement->getMaterial(), Event::CELL, _ptype, 0, 1, _state);
			_pElement->_ABalance--;
		}
		else
		{
			_pElement->decBAtoms();
                        pSub->_evLog.performed(_pElement->getMaterial(), Event::CELL, _ptype, 0, 3, _state);
			_pElement->_BBalance--;
		}
	}
	if (Domains::global()->PM()->getMaterial(_pElement->getMaterial())._selfDiffusion && 
            Domains::global()->PM()->getIorV(_pElement->getMaterial(), _ptype) == Kernel::IV_V)
	{
		if(_pDomain->_rng_dom.rand() > _pElement->getAlloyFraction())
		{
			_pElement->incAAtoms();
                        pSub->_evLog.performed(_pElement->getMaterial(), Event::CELL, _ptype, 0, 0, _state);
			_pElement->_ABalance++;
		}
		else
		{
			_pElement->incBAtoms();
                        pSub->_evLog.performed(_pElement->getMaterial(), Event::CELL, _ptype, 0, 2, _state);
			_pElement->_BBalance++;
		}
	}
	_pDomain->_pRM->remove(this, _pElement);
	_pDomain->_pMesh->remove(this);

	if (Domains::global()->IsPoisson())
		_pDomain->_pMesh->getParticleToNodeHandler()->remove(this);
}

/* Takes the axis of the defect.
 * Produces a permutation on it, and also changes the sign
 * Rotates it according to the current wafer rotation.
 */
void MobileParticle::setAxes(Kernel::SubDomain *pSub)
{
	//changes the sign
	for(int j=0; j<3; ++j) //x, y or z
	{
		int sign = (pSub->_rng.rand() < 0.5? 1:-1);
		_axes[j] = sign * _pDomain->_pMPPar->_axes[_pElement->getMaterial()][_ptype][_state][j];
	}
	//permutation. There are six possible ones abc acb bac bca cab cba
	int permutations = int(pSub->_rng.rand()*6.)+1;
	vector<unsigned> perm;
	perm.push_back(0);	perm.push_back(1); perm.push_back(2);
	for(int i=1; i<permutations; ++i)
		std::next_permutation(perm.begin(), perm.end());
	//perm contains now the permuted indexes.
	Coordinates temp; //otherwise it can rewrite itself.
	for(int j=0; j<3; ++j)
		temp[j] = _axes[perm[j]];
	_axes = temp;
	if(_axes.abs())
	{
		_axes = _pDomain->_pLat[_pElement->getMaterial()]->toRotated(_axes);
		//normalize axes
		_axes /= _axes.abs();
	}
}

float MobileParticle::getRate(unsigned ev, float kT) const
{
	Kernel::M_TYPE mt = _pElement->getMaterial();
#ifdef NUMODEL
	unsigned migEv = 5;
	if(ev < 6) //migration
	{
		Kernel::MeshElement *pEle = _pElement;
		_longHopFactor =  _pDomain->_pMesh->longHopFactor(_pElement->getIndex());
		Coordinates orig_delta(0,0,0);
		const float l = _pDomain->_pMePar->_lambda[pEle->getMaterial()];
		if(ev == 0)
			orig_delta._x = l;
		else if(ev == 1)
			orig_delta._x = -l;
		else if(ev == 2)
			orig_delta._y = l;
		else if(ev == 3)
			orig_delta._y = -l;
		else if(ev == 4)
			orig_delta._z = l;
		else
			orig_delta._z = -l;
		Coordinates to = _coord + (orig_delta * _longHopFactor);
		Kernel::M_TYPE mt = _pElement->getMaterial();
		Mesh::JUMP_ACTIONS jmp = _pDomain->_pMesh->checkMove(pEle, to);
		IO::Arrhenius arr = _pDomain->_pMPPar->_arr[mt][_ptype][_state][ev](_pElement);
		double x1 = _pElement->getAlloyFraction();
		double x2 = pEle->getAlloyFraction();
		char iorv = Domains::global()->PM()->getIorV(mt, _ptype);

		if(iorv == Kernel::IV_I)
			_AorB[ev] = _pDomain->_rng_dom.rand() > _pElement->getAlloyFraction() ? true : false;
		else if(iorv == Kernel::IV_V)
			_AorB[ev] = _pDomain->_rng_dom.rand() > pEle->getAlloyFraction() ? true : false;
		else
			ERRORMSG("MobileParticle not I or V is trying to get a rate for moving lattice atoms");
		if(jmp == Mesh::JUMP_OK && x1 != x2)
		{
			// Workaround for unresolved singularities [ln at 0]
			if(x1 == 0)
			{
				if(iorv == Kernel::IV_I && !_AorB[ev])
					return 0;
				x1 = 0.001;
			}
			if(x2 == 0)
			{
				if(iorv == Kernel::IV_V && !_AorB[ev])
					return 0;
				x2 = 0.001;
			}
			if(x1 == 1)
			{
				if(iorv == Kernel::IV_I && _AorB[ev])
					return 0;
				x1 = 0.999;
			}
			if(x2 == 1)
			{
				if(iorv == Kernel::IV_V && _AorB[ev])
					return 0;
				x2 = 0.999;
			}
			if(_AorB[ev])
				arr._ener -= .5 * (_pDomain->_pAlPar->nu(iorv==Kernel::IV_I?x2:x1, kT) - _pDomain->_pAlPar->nu(iorv==Kernel::IV_I?x1:x2, kT));
			else
				arr._ener += .5 * (_pDomain->_pAlPar->nu(iorv==Kernel::IV_I?x2:x1, kT) - _pDomain->_pAlPar->nu(iorv==Kernel::IV_I?x1:x2, kT));
		}
		double rate = arr.getRate(kT);
		return rate / (_longHopFactor * _longHopFactor);
	}
#else
	unsigned migEv = 0;
	if(ev == 0) //migration
	{
		double rate = _pDomain->_pMPPar->_arr[mt][_ptype][_state][ev](_pElement).getRate(kT);
		_longHopFactor =  _pDomain->_pMesh->longHopFactor(_pElement->getIndex());
		return rate / (_longHopFactor * _longHopFactor);
	}
#endif
	else if(ev == (migEv + 1) || ev == (migEv + 2)) //breakup prefactors not considered yet in stateChanges
	{
		double stateFormation = _pDomain->_pSM->bindingShift(_ptype, (ev == 1? POS_0 : POS_1), _state, _pElement);
		if(_pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._isCluster)
			return breakUpClusterRate(kT, ev);
		return (_pDomain->_pMPPar->_arr[mt][_ptype][_state][ev  ](_pElement).getRate(kT) +
				_pDomain->_pMPPar->_arr[mt][_ptype][_state][ev+5](_pElement).getRate(kT))
				* std::exp(-stateFormation/kT);
	}
	else if(ev == (migEv + 3) || ev == (migEv + 4)) //FT break-up
	{
		return _pDomain->_pMPPar->_arr[mt][_ptype][_state][ev  ](_pElement).getRate(kT) +
				_pDomain->_pMPPar->_arr[mt][_ptype][_state][ev+5](_pElement).getRate(kT);
	}
	else if (ev == (migEv + 5)) //update
	{
		return _pDomain->_pMPPar->_arr[mt][_ptype][0][ev](_pElement).getRate(kT);
	}
	else
		ERRORMSG("Event " << ev << " not defined in MobileParticle");
	return 0;
}

//computes the rate when the break up particles creates a cluster (i.e., VSi -> VC^CSi) (C -> AB )
// The rate is computed as clToPartRate / clToSplitRate * splitRate = AB->C * C->A+B  / AB -> A+B
// clToPartRate  = rate Cluster AB to recombine into part C
// splitRate     = rate particle C to break up into particles A + B
// clToSplitRate = rate Cluster AB to break into A + B
double MobileParticle::breakUpClusterRate(float kT, int ev) const
{
	Kernel::M_TYPE mt = _pElement->getMaterial();
	double splitRate = _pDomain->_pMPPar->_arr[mt][_ptype][_state][ev  ](_pElement).getRate(kT) +
			_pDomain->_pMPPar->_arr[mt][_ptype][_state][ev+5](_pElement).getRate(kT);
	P_TYPE emit0 = _pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[0];
	P_TYPE emit1 = _pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[1];
	const map<unsigned, float> &types = _pDomain->_pMPPar->_mctype[mt][emit0][emit1];
	if(types.size() != 1)
		ERRORMSG("Cannot break up a particle into more than one type of cluster");
	unsigned edtype = types.begin()->first;
	const EDType * clParams = _pDomain->_pClPar->getParams(mt, edtype);
	ID thisMap = ClusterParam::pt2ID(mt, emit0);
	Cluster::addMap(thisMap, ClusterParam::pt2ID(mt, emit1), clParams->_ivmodel);
	const EDType::CLType *myParams = _pDomain->_pClPar->getParams(edtype, thisMap);

	double clToPartRate  = Cluster::recombinationRate(getElement(), edtype, myParams, thisMap, kT);
	double clToSplitRate =
			Cluster::emissionRate(getElement(), edtype, myParams, thisMap, emit0, kT) +
			Cluster::emissionRate(getElement(), edtype, myParams, thisMap, emit1, kT);
	if(clToSplitRate)
		return clToPartRate / clToSplitRate * splitRate;
	else
		return 0;
}

//move particle, look for neighbors, interact with incoming particles.
void MobileParticle::perform (Kernel::SubDomain *pSub, unsigned eventType)
{
	if(eventType != 0)
		pSub->_evLog.performed(_pElement->getMaterial(), Event::MOBILEPARTICLE, _ptype, 0, eventType, _state);
	switch(eventType)
	{
#ifdef NUMODEL
	case 0:
		return migrate(pSub, 0);
	case 1:
		return migrate(pSub, 1);
	case 2:
		return migrate(pSub, 2);
	case 3:
		return migrate(pSub, 3);
	case 4:
		return migrate(pSub, 4);
	case 5:
		return migrate(pSub, 5);
	case 6:
		return breakup(pSub, POS_0);
	case 7:
		return breakup(pSub, POS_1);
	case 8:
		return emit(pSub, Kernel::IV_I);
	case 9:
		return emit(pSub, Kernel::IV_V);
	case 10:
		return updateState(pSub);
#else
	case 0:
		return migrate(pSub);
	case 1:
		return breakup(pSub, POS_0);
	case 2:
		return breakup(pSub, POS_1);
	case 3:
		return emit(pSub, Kernel::IV_I);
	case 4:
		return emit(pSub, Kernel::IV_V);
	case 5:
		return updateState(pSub);
#endif
	default:
		ERRORMSG("Unknown case in MobileParticle::perform");
		break;
	}
}

void MobileParticle::emit(Kernel::SubDomain *pSub, unsigned char iorv)
{
	IO::ParameterManager *pPM = Domains::global()->PM();
	P_POS pos = pPM->getPPos(_ptype);
	Kernel::M_TYPE mt = _pElement->getMaterial();
	const bool binary = pPM->getMaterial(mt)._binary;
	P_POS otherPos = (pos == POS_0 && binary)? POS_1 : POS_0;
	P_TYPE fam = pPM->getFamily(_ptype);
	if(pos == Kernel::POS_I || pos == Kernel::POS_V || fam == Kernel::V_TYPE) //only substitutionals
		ERRORMSG("Incorrect emission in mobile particle");
	P_TYPE emit1 = pPM->getParticle(mt, fam,  (iorv == Kernel::IV_I? Kernel::POS_V:Kernel::POS_I));
	P_TYPE emit2 = pPM->iorv2pt    (mt, iorv, (iorv == Kernel::IV_I? otherPos: pos));
	assert(emit1 != Kernel::UNDEFINED_TYPE && emit2 != Kernel::UNDEFINED_TYPE);

	Kernel::Coordinates newCoords(_coord);
	breakPosition(pSub, newCoords, _pElement, _pDomain->_pMPPar->_interact_radius[mt][emit1][emit2]);
	if(Domains::global()->PM()->getIorV(mt, _ptype) == Kernel::IV_NONE) //new type is still a mobilepart
		new MobileParticle(pSub, emit2, _state, _pDomain, newCoords, _orig); //self inserts
	else //implement Cluster
		ERRORMSG("Cluster creation not implemented in MobileParticle yet");

	_pDomain->_pRM->remove(this, _pElement);
	_ptype = emit1;
	_pDomain->_pRM->insert(this, _pElement);
	if (Domains::global()->IsPoisson())
		_pDomain->_pMesh->getParticleToNodeHandler()->insert(this);
}

#ifdef NUMODEL
void MobileParticle::migrate(Kernel::SubDomain *pSub, unsigned ev)
{
	Kernel::MeshElement *pEle = _pElement;
	Coordinates oldc = _coord;
	Coordinates oldorig = _orig;
	Coordinates orig_delta(0,0,0);
	const float l = _pDomain->_pMePar->_lambda[pEle->getMaterial()];
	if(_axes.abs() == 0)
	{
		if(ev == 0)
			orig_delta._x = l;
		else if(ev == 1)
			orig_delta._x = -l;
		else if(ev == 2)
			orig_delta._y = l;
		else if(ev == 3)
			orig_delta._y = -l;
		else if(ev == 4)
			orig_delta._z = l;
		else
			orig_delta._z = -l;
	}
	else
		orig_delta = _axes * (pSub->_rng.rand() < 0.5? l:-l);
	const Coordinates delta(orig_delta*_longHopFactor);

	std::pair<P_TYPE, unsigned> stateID = std::make_pair(_ptype, _state);   //correction to make state structure fit into jump actions.
	Kernel::M_TYPE mt = _pElement->getMaterial();
	P_TYPE pt = _ptype;
	unsigned st = _state;
	unsigned uLHF = _longHopFactor;
	Mesh::JUMP_ACTIONS jmp = _pDomain->_pMesh->jumpPosition(_coord, delta, pEle, &stateID, _orig, _longHopFactor*_longHopFactor);
	vector<Particle *> parts;
	if (jmp != Mesh::JUMP_OK)
	{
		pSub->_evLog.performed(mt, Event::MOBILEPARTICLE, _ptype, 0, 12, _state);
		return;
	}

	if(_pElement != pEle ) //change of elements
	{
		if(Domains::global()->PM()->isAmorphous(pEle->getMaterial()))
		{
			if(!Domains::global()->PM()->isImpurity(mt,_ptype))
			{
				HIGHMSG("Deleting " << Domains::global()->PM()->getParticleName(mt,_ptype));
				delete this;
			}
			else
			{
				HIGHMSG("Desorption in a/c interface of " << Domains::global()->PM()->getParticleName(mt,_ptype) << ". Removing "
						<< Domains::global()->PM()->getParticleName(mt,_pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[1])
						<< " and leaving " << Domains::global()->PM()->getParticleName(mt,_pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[0] ));
				P_TYPE emit  = _pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[0];
				_pDomain->_pRM->remove(this, _pElement);
				_ptype = emit;
				_state = 0;
				_pDomain->_pRM->insert(this, _pElement);
			}
			return;
		}
		INTERFACE_ACTIONS result = interactSurface(pSub, _pElement, pEle);
		if(result == INTERACTION_REJECTED)
			goto rejected;
		if(result == INTERACTION_EXECUTED)
			goto executed;
		if ( !_pDomain->_pSM->checkStateBarrier(pSub, _ptype, _state , _pElement ,pEle))
			goto rejected;

		double newLHF = _pDomain->_pMesh->longHopFactor(pEle->getIndex());  // Martin-Bragado et al. SSE . 155-155 (2008) 202-206
		if(_longHopFactor < newLHF && pSub->_rng.rand() > _longHopFactor/double(newLHF))  //long - short hop probability
			goto rejected;

		if(newLHF != _longHopFactor)
		{
			unsigned ix, iy, iz;
			_pDomain->_pMesh->getIndicesFromIndex(pEle->getIndex(), ix, iy, iz);

			if(delta._x > 0)
				_coord._x = _pDomain->_pMesh->getLines(0)[ix]   + pSub->_rng.rand()*newLHF*orig_delta._x;
			else if(delta._x < 0)
				_coord._x = _pDomain->_pMesh->getLines(0)[ix+1] + pSub->_rng.rand()*newLHF*orig_delta._x;
			if(delta._y > 0)
				_coord._y = _pDomain->_pMesh->getLines(1)[iy]   + pSub->_rng.rand()*newLHF*orig_delta._y;
			else if(delta._y < 0)
				_coord._y = _pDomain->_pMesh->getLines(1)[iy+1] + pSub->_rng.rand()*newLHF*orig_delta._y;
		}

		// Self-Diffusion
		char iorv = Domains::global()->PM()->getIorV(mt, _ptype);
		if(iorv == Kernel::IV_I)
		{
			if(_AorB[ev])
			{
				_pElement->decAAtoms();
				pEle->incAAtoms();
			}
			else
			{
				_pElement->decBAtoms();
				pEle->incBAtoms();
			}
		}
		else if (iorv == Kernel::IV_V)
		{
			if(_AorB[ev])
			{
				_pElement->incAAtoms();
				pEle->decAAtoms();
			}
			else
			{
				_pElement->incBAtoms();
				pEle->decBAtoms();
			}
		}
		else
			ERRORMSG("Only interstitials and vacancies are allowed to move lattice atoms");
		assert(Domains::global()->PM()->isParticleDefined(_ptype, mt));
	}
	_pDomain->_pMesh->remove(this);
	_pDomain->_pMesh->insert(this, pEle);
	_pDomain->_pRM->remove(this, _pElement);
	_pDomain->_pRM->insert(this, pEle);
	pSub->_evLog.performed(mt, Event::MOBILEPARTICLE, pt, 0, (uLHF == 1? ev : 11), st);
#else
void MobileParticle::migrate(Kernel::SubDomain *pSub)
{
	Kernel::MeshElement *pEle = _pElement;
	Coordinates oldc = _coord;
	Coordinates oldorig = _orig;
	Coordinates orig_delta(0,0,0);
	const float l = _pDomain->_pMePar->_lambda[pEle->getMaterial()];
	if(_axes.abs() == 0)
	{
		const int axis = int(pSub->_rng.rand() * 6.);
		if(axis == 0)
			orig_delta._x = l;
		else if(axis == 1)
			orig_delta._x = -l;
		else if(axis == 2)
			orig_delta._y = l;
		else if(axis == 3)
			orig_delta._y = -l;
		else if(axis == 4)
			orig_delta._z = l;
		else
			orig_delta._z = -l;
	}
	else
		orig_delta = _axes * (pSub->_rng.rand() < 0.5? l:-l);
	const Coordinates delta(orig_delta*_longHopFactor);

	std::pair<P_TYPE, unsigned> stateID = std::make_pair(_ptype, _state);   //correction to make state structure fit into jump actions.
	Kernel::M_TYPE mt = _pElement->getMaterial();
	P_TYPE pt = _ptype;
	unsigned st = _state;
	unsigned uLHF = _longHopFactor;
	Mesh::JUMP_ACTIONS jmp = _pDomain->_pMesh->jumpPosition(_coord, delta, pEle, &stateID, _orig, _longHopFactor*_longHopFactor);
	vector<Particle *> parts;
	if (jmp != Mesh::JUMP_OK)
	{
		pSub->_evLog.performed(mt, Event::MOBILEPARTICLE, _ptype, 0, 7, _state);
		return;
	}

	if(_pElement != pEle ) //change of elements
	{
		if(Domains::global()->PM()->isAmorphous(pEle->getMaterial()))
		{
			if(!Domains::global()->PM()->isImpurity(mt,_ptype))
			{
				HIGHMSG("Deleting " << Domains::global()->PM()->getParticleName(mt,_ptype));
				delete this;
			}
			else
			{
				HIGHMSG("Desorption in a/c interface of " << Domains::global()->PM()->getParticleName(mt,_ptype) << ". Removing "
						<< Domains::global()->PM()->getParticleName(mt,_pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[1])
						<< " and leaving " << Domains::global()->PM()->getParticleName(mt,_pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[0] ));
				P_TYPE emit  = _pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[0];
				_pDomain->_pRM->remove(this, _pElement);
				_ptype = emit;
				_state = 0;
				_pDomain->_pRM->insert(this, _pElement);
			}
			return;
		}
		INTERFACE_ACTIONS result = interactSurface(pSub, _pElement, pEle);
		if(result == INTERACTION_REJECTED)
			goto rejected;
		if(result == INTERACTION_EXECUTED)
			goto executed;
		if ( !_pDomain->_pSM->checkStateBarrier(pSub, _ptype, _state , _pElement ,pEle))
			goto rejected;

		double newLHF = _pDomain->_pMesh->longHopFactor(pEle->getIndex());  // Martin-Bragado et al. SSE . 155-155 (2008) 202-206
		if(_longHopFactor < newLHF && pSub->_rng.rand() > _longHopFactor/double(newLHF))  //long - short hop probability
			goto rejected;

		if(newLHF != _longHopFactor)
		{
			unsigned ix, iy, iz;
			_pDomain->_pMesh->getIndicesFromIndex(pEle->getIndex(), ix, iy, iz);

			if(delta._x > 0)
				_coord._x = _pDomain->_pMesh->getLines(0)[ix]   + pSub->_rng.rand()*newLHF*orig_delta._x;
			else if(delta._x < 0)
				_coord._x = _pDomain->_pMesh->getLines(0)[ix+1] + pSub->_rng.rand()*newLHF*orig_delta._x;
			if(delta._y > 0)
				_coord._y = _pDomain->_pMesh->getLines(1)[iy]   + pSub->_rng.rand()*newLHF*orig_delta._y;
			else if(delta._y < 0)
				_coord._y = _pDomain->_pMesh->getLines(1)[iy+1] + pSub->_rng.rand()*newLHF*orig_delta._y;
		}

		if(_pElement->getBAtoms() != 0 || pEle->getBAtoms() != 0)
		{
			float kT = _pDomain->_pRM->getkT();
			float T = _pDomain->_pRM->getT();
			// Segregation and concentrations
			IO::Arrhenius migFrom = _pDomain->_pMPPar->_arr[mt][_ptype][_state][0](_pElement);
			IO::Arrhenius formFrom = _pDomain->_pMPPar->_form[mt][_ptype](_pElement);
			IO::Arrhenius migTo = _pDomain->_pMPPar->_arr[mt][_ptype][_state][0](pEle);
			IO::Arrhenius formTo = _pDomain->_pMPPar->_form[mt][_ptype](pEle);
			float x1 = _pElement->getEffectiveAlloyFraction();
			float x2 = pEle->getEffectiveAlloyFraction();

			float sign = Domains::global()->PM()->getIorV(mt, _ptype) ? -1 : 1;
			float effFormEnerFrom = formFrom._ener + sign * (_pDomain->_pAlPar->getMixingEnergy(mt, x1, T) +
					(.5 - x1) * _pDomain->_pAlPar->getDerMixingEnergy(mt, x1, T));
			float effFormEnerTo = formTo._ener + sign * (_pDomain->_pAlPar->getMixingEnergy(mt, x2, T) +
					(.5 - x2) * _pDomain->_pAlPar->getDerMixingEnergy(mt, x2, T));
			double QFrom = effFormEnerFrom + migFrom._ener;
			double QTo = effFormEnerTo + migTo._ener;

			double pBarrier = (formTo._pref / formFrom._pref) * (migTo._pref / migFrom._pref) * std::exp( -(QTo - QFrom) / kT);
			if (pBarrier < 1)
				if (pSub->_rng.rand() > pBarrier)
					goto rejected;
		}

		//Self-Diffusion
		selfdiffusion(pSub, pEle);

		assert(Domains::global()->PM()->isParticleDefined(_ptype, mt));
		_pDomain->_pRM->  remove(this, _pElement);
		_pDomain->_pMesh->remove(this);
		_pDomain->_pMesh->insert(this, pEle);
		_pDomain->_pRM->  insert(this, pEle);
	}
	pSub->_evLog.performed(mt, Event::MOBILEPARTICLE, pt, 0, (uLHF == 1? 0 : 6), st);
#endif

	//order is very important, because the scheduler uses the element to pick up the subdomain.
	if(uLHF == 1)
		_pDomain->_pMesh->getInteractions(pSub, this, parts);
	if(parts.size())
	{
		Particle *intPar = parts[int(pSub->_rng.rand() * parts.size())];
		std::vector<Particle *> dummy;
		intPar->getDefect()->interact(pSub, this, intPar, dummy);
	}
	return;

rejected:
	_orig = oldorig;
	_coord = oldc;
#ifdef NUMODEL
	pSub->_evLog.performed(mt, Event::MOBILEPARTICLE, pt, 0, 12, st);
#else
	pSub->_evLog.performed(mt, Event::MOBILEPARTICLE, pt, 0, 7, st);
#endif
	return;

executed:
#ifdef NUMODEL
	pSub->_evLog.performed(mt, Event::MOBILEPARTICLE, pt, 0, (uLHF == 1? ev : 11), st);
#else
	pSub->_evLog.performed(mt, Event::MOBILEPARTICLE, pt, 0, (uLHF == 1? 0 : 6), st);
#endif
	return;
}

void MobileParticle::breakup(Kernel::SubDomain *pSub, P_POS pos)
{
	Kernel::M_TYPE mt = _pElement->getMaterial();
	const bool toCluster = _pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._isCluster;
	P_TYPE emit  = _pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[0];
	P_TYPE emit2 = _pDomain->_pMPPar->_breakUp[mt][_ptype][_state]._emit[1];
	unsigned stateIorV = _pDomain->_pSM->breakUpState(_ptype, pos, _state, _pElement->getMaterial());
	_pDomain->_pRM->remove(this, _pElement);
	_ptype = emit;
	_state = 0;
	if (Domains::global()->IsPoisson())
		_pDomain->_pMesh->getParticleToNodeHandler()->insert(this);
	_pDomain->_pRM->insert(this, _pElement);

	Kernel::Coordinates newCoords(_coord);
	breakPosition(pSub, newCoords, _pElement, _pDomain->_pMPPar->_interact_radius[mt][emit][emit2]);
	MobileParticle * pMP = new MobileParticle(pSub, emit2, stateIorV, _pDomain, newCoords, newCoords); //self inserts
	pMP->updateState(pSub);
	if(toCluster)
	{
		std::vector<Particle *> dummy;
		internalInteract(pSub, pMP, true, dummy); //forces the creation of a cluster
	}
}

//implements the particle interaction
//in principle, interactions should not be discarded here. Use canInteract for that.
Defect * MobileParticle::interact(Kernel::SubDomain *pSub, Defect *dyn, Particle *, std::vector<Particle *> &parts)
{
	return internalInteract(pSub, dyn, false, parts);
}

Defect * MobileParticle::internalInteract(Kernel::SubDomain *pSub, Defect *def, bool bForce, std::vector<Particle *> &parts)
{
	Kernel::M_TYPE mt = _pElement->getMaterial();
	switch(def->getEType())
	{
	case MOBILEPARTICLE:
	{
		MobileParticle *mp = dynamic_cast<MobileParticle *>(def);
		switch(_pDomain->_pMPPar->_interact[mt][int(_ptype)][int(mp->getPType())]) //result
				{
		case MOBILEPARTICLE:
		{
			P_TYPE newPt = _pDomain->_pMPPar->_int_result[int(_ptype)][mp->getPType()];
			unsigned resState = _pDomain->_pSM->interactionState(_ptype, _state, mp->getPType(), mp->getState(), _pElement, newPt);
			if(resState == Kernel::UNDEFINED_STATE)
				return this; //no reaction is performed if charge is out of bounds.
			{  //check energetics
				float initEner = _pDomain->_pMPPar->_form[mt][_ptype](_pElement)._ener +
						_pDomain->_pMPPar->_form[mt][mp->getPType()](_pElement)._ener;
				float endEner = _pDomain->_pMPPar->_form[mt][newPt](_pElement)._ener;
				float diffEner = endEner - initEner;
				if(bForce == false && diffEner >0 && pSub->_rng.rand() > exp(-diffEner/_pDomain->_pRM->getkT()))
					return this;
			}
			//coordinates are the incoming ones, to avoid Moire effects
			pSub->_reLog.reaction(_pElement->getMaterial(), MOBILEPARTICLE, _ptype,_state, MOBILEPARTICLE, mp->getPType(), mp->getState());
			if(Domains::global()->PM()->getIorV(mt, _ptype) == Kernel::IV_NONE) //if I am the impurity
			{
				_pDomain->_pRM->remove(this, _pElement);
				_ptype = newPt;
				_state = resState;
				if (Domains::global()->IsPoisson())
					_pDomain->_pMesh->getParticleToNodeHandler()->insert(this);
				_pDomain->_pRM->insert(this, _pElement);
				delete mp;
				return this;
			}
			else  //mp is the impurity
			{
				_pDomain->_pRM->remove(mp, mp->_pElement);
				mp->_ptype = newPt;
				mp->_state = resState;
				if (Domains::global()->IsPoisson())
					_pDomain->_pMesh->getParticleToNodeHandler()->insert(mp);
				_pDomain->_pRM->insert(mp, mp->_pElement);
				delete this;
				return mp;
			}
		}
		case CLUSTER:
		{
			double rn = pSub->_rng.rand();
			unsigned edtype=0;
			const map<unsigned, float> &types = _pDomain->_pMPPar->_mctype[_pElement->getMaterial()][int(_ptype)][mp->getPType()];
			//compute probability for correct extended defect
			for(map<unsigned, float>::const_iterator it=types.begin(); it!=types.end(); ++it)
			{
				if(rn < it->second)
				{
					edtype = it->first;
					rn -= it->second;
					break;
				}
				else
					rn -= it->second;
			}
			if(rn > 0)
			{
				for(map<unsigned, float>::const_iterator it=types.begin(); it!=types.end(); ++it)
					LOWMSG(_pDomain->_pClPar->getParams(_pElement->getMaterial(),it->first)->_name << " " << it->second);
				ERRORMSG("MobileParticle: Sorry, but could not find a proper cluster for interaction. " << rn);
			}
			if(bForce || _pDomain->_pClPar->reactionPossible(pSub, _pDomain, _pElement, edtype, mp->getPType(), _ptype, _pDomain->_pRM->getkT()))
			{
				pSub->_reLog.reaction(_pElement->getMaterial(), MOBILEPARTICLE, _ptype,_state, MOBILEPARTICLE, mp->getPType(), mp->getState());
				Kernel::Domain *pD = _pDomain; //because next lines deletes 'this' pointer
				Cluster *pCL = new Cluster(pSub, _pDomain, edtype, mp, this); //cluster removes the mobile particles
				pD->_pRM->insert(pCL, pCL->getElement());
				return pCL;
			}
			return this;
		}
		case EMPTY:
			pSub->_reLog.reaction(_pElement->getMaterial(), MOBILEPARTICLE, _ptype,_state, MOBILEPARTICLE, mp->getPType(), mp->getState());
			delete mp;
			delete this;
			return 0;
		case SINK:
			pSub->_reLog.reaction(_pElement->getMaterial(), MOBILEPARTICLE, _ptype,_state, MOBILEPARTICLE, mp->getPType(), mp->getState());
			delete mp;
			return 0;
		default:
			ERRORMSG("MobileParticle: Event interaction " <<
					Event::getEName(Event::E_TYPE(_pDomain->_pMPPar->_interact[_pElement->getMaterial()][int(_ptype)][int(mp->getPType())])) <<
					" not supported yet.");
				}
		return this;
	}
	default:
		return def->interact(pSub, this, 0, parts);
	}
}

//this canInteract does not check the energetics. That is done in the real interaction.
bool MobileParticle::canInteract(Kernel::SubDomain *pSub, const Particle *dyn, const Particle *sta) const
{
	Defect *def = dyn->getDefect();
	assert(dyn->getDefect() != this);
	assert(sta != dyn);
	if(Domains::global()->PM()->isAmorphous(_pElement->getMaterial()) || Domains::global()->PM()->isAmorphous(def->getElement()->getMaterial()))
		return false;
	switch(def->getEType())
	{
	case MOBILEPARTICLE:
	{
		Kernel::M_TYPE mt = _pElement->getMaterial();
		const MobileParticle *mp = static_cast<const MobileParticle *>(def);
		if(_pDomain->_pMPPar->_canInteract[mt][int(_ptype)][int(mp->getPType())]  == false ||
				_pDomain->_pSM->canInteract(_ptype, _state, mp->_ptype, mp->_state, _pElement) == false )
			return false;
		// compute neighbor distance
		Coordinates n = dyn->getCoordinates();
		_pDomain->_pMesh->setPeriodicRelative(_coord, n);
		const float cdist2 = n*n;
		return (_pDomain->_pMPPar->_interact_radiusSq[mt][_ptype][mp->_ptype] > cdist2);
	}
	default:
		return def->canInteract(pSub, sta, dyn);
	}
}

void MobileParticle::updateState(Kernel::SubDomain *pSub)
{
	unsigned state = _pDomain->_pSM->changeState(pSub, _ptype,_pElement);
	if(state != _state)
	{
		_pDomain->_pRM->remove(this, _pElement);
		_state = state;
		_pDomain->_pRM->insert(this, _pElement);
	}
}

ID MobileParticle::getID() const
{
	return ClusterParam::pt2ID(_pElement->getMaterial(), _ptype);
}

void MobileParticle::deletePart(Kernel::SubDomain *, std::vector<Particle *> &parts)
{
	assert(parts.back() == this);
	parts.pop_back();
	delete this;
}

}


