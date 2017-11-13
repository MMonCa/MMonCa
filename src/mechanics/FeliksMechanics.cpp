/*
 * FELIKSMechanics.cpp
 *
 *  Created on: May 8, 2013
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

#include "FeliksMechanics.h"

#ifdef MMONCA_FELIKS
#include "io/MaterialProperties.h"
#include "io/ParameterManager.h"
#include "kernel/Domain.h"
#include "kernel/RateManager.h"
#include "kernel/Mesh.h"
#include "kernel/MeshElement.h"

#include "okmc/Defect.h"

#include "feliks/Analysis/FEanalysis/festatic.h"
#include "feliks/Io/io.h"
#include "feliks/Materials/material.h"

// global definitions
double          global_macheps;
bool            global_bigEndianSystem;

using std::vector;

extern "C" {double dlamch_(char *, short );}

namespace Mechanics {

static void initialize();

FELIKSMechanics::FELIKSMechanics(Kernel::Domain *p) : MechInterface(p)
{
	// open the files where the programs logs the information of the analysis
	LOWMSG2("Initializing FELIKS.");
	logger      :: openLogFiles();
	factorcombo :: initializePropfactorCombinations();
	initialize();

   // vertices are (nx+1)(ny+1)(nz+1)
	const vector<float> &lx = _pDomain->_pMesh->getLines(0);
	vector<float> ly;
	vector<float> lz;
	if(_dim > 1)
		ly = _pDomain->_pMesh->getLines(1);
	else
	{
		ly.push_back(0);
		ly.push_back(1);
	}

	if(_dim > 2)
		lz = _pDomain->_pMesh->getLines(2);
	else
	{
		lz.push_back(0);
		lz.push_back(1);
	}
	_nX = lx.size();
	_nY = ly.size();
	_nZ = lz.size();

	_vertices.reserve( _nX*_nY*_nZ*3);
	for(unsigned i=0; i<_nX; ++i)
		for(unsigned j=0; j<_nY; ++j)
			for(unsigned k=0; k<_nZ; ++k)
			{
				_vertices.push_back(lx[i]);
				_vertices.push_back(ly[j]);
				_vertices.push_back(lz[k]);
			}
	// elements _connectivity
	_connectivity.reserve((_nX-1)*(_nY-1)*(_nZ-1)*8);
	for(unsigned i=0; i<_nX-1; ++i)
			for(unsigned j=0; j<_nY-1; ++j)
				for(unsigned k=0; k<_nZ-1; ++k)
				{
					_connectivity.push_back(idx2vert(i  ,j  ,k  ));
					_connectivity.push_back(idx2vert(i+1,j  ,k  ));
					_connectivity.push_back(idx2vert(i+1,j+1,k  ));
					_connectivity.push_back(idx2vert(i  ,j+1,k  ));
					_connectivity.push_back(idx2vert(i  ,j  ,k+1));
					_connectivity.push_back(idx2vert(i+1,j  ,k+1));
					_connectivity.push_back(idx2vert(i+1,j+1,k+1));
					_connectivity.push_back(idx2vert(i  ,j+1,k+1));
				}
	// boundary _conditions
	for(unsigned i=0; i<_nX; ++i)
			for(unsigned j=0; j<_nY; ++j)
				for(unsigned k=0; k<_nZ; ++k)
					if(i == 0 || i == _nX -1)
						_bc[idx2vert(i,j,k)] = 1;
					else if (j == 0 || j == _nY -1)
						_bc[idx2vert(i,j,k)] = 2;
					else if (k == 0 || k == _nZ -1)
						_bc[idx2vert(i,j,k)] = 4;
	// materials
	IO::ParameterManager *pPM = Domains::global()->PM();
	for(Kernel::M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
	{
		std::string name = pPM->getMaterialName(mt);
		const IO::MaterialProperties &m = pPM->getMaterial(mt);
		_material2idx[mt] = _materials.size();
		_materials.push_back( feliks::mmoncaMaterial(name, m.thermalExp_, m.tempRef_, m.E_, m.nu_,
				m.eigenstrains_) );
		if(m._alloyName != "none") //there is an alloy
		{
			name = m._alloyName;
			_materials.push_back( feliks::mmoncaMaterial(name, m.alloythermalExp_, m.alloytempRef_, m.alloyE_, m.alloynu_,
					m.eigenstrains_) );
		}
	}

	for(unsigned nE = 0; nE < _pDomain->_pMesh->size(); ++nE)
	{
		Kernel::MeshElement * pME = _pDomain->_pMesh->getElement(nE);
		unsigned ii, jj, kk;
		_pDomain->_pMesh->getIndicesFromIndex(nE, ii, jj, kk);
		if(_dim < 2 && jj != 0)
			continue;
		if(_dim < 3 && kk != 0)
			continue;
		// element state
		// for element 1: composition + state
		Kernel::M_TYPE mt = pME->getMaterial();
		vector<double> comp_el(_materials.size(), 0.);
		const IO::MaterialProperties &element = pPM->getMaterial(mt);
		float comp = pME->getBAtoms();
		if(comp != 0)
			comp = comp / pME->getVolume() / element._densityCm3 * 1e21; // atm/(nm)^3
		comp_el[_material2idx[mt]] = 1 - comp;
		if(comp != 0)
			comp_el[_material2idx[mt] + 1] = comp;
		feliks::elementState es(_materials.size());
		es.setState(comp_el, _pDomain->_pRM->getT(), 0.0, 0.0, 0.0, map<string, unsigned>());     //  temp, eps1, eps2, eps3
		_elstates.push_back(es);
	}
	// Create the interface. Do not do any computation
	_pmi = new feliks::mmoncaInterface(_vertices, _connectivity, _elstates, _bc, _materials);
	_theAnalysis = new staticFEanalysis(*_pmi);
	LOWMSG("..Done.");
}

//the vertex with i,j,k is at i(YZ) + jZ + k
unsigned FELIKSMechanics::idx2vert(unsigned i, unsigned j, unsigned k) const
{
	return i*_nY*_nZ + j*_nZ + k;
}

FELIKSMechanics::~FELIKSMechanics()
{
	delete _theAnalysis;
	delete _pmi;
	// Deallocate global lists
	logger      :: closeLogFiles();
	factorcombo :: freePropfactorCombinations();
	material    :: freeMaterials();
	eltype      :: freeEltypes();
}

void FELIKSMechanics::import() const
{
	LOWMSG2("Computing stress and strain using FELIKS...");
	unsigned i=0;
	for(unsigned nE = 0; nE < _pDomain->_pMesh->size(); ++nE)
	{
		Kernel::MeshElement * pME = _pDomain->_pMesh->getElement(nE);
		unsigned ii, jj, kk;
		_pDomain->_pMesh->getIndicesFromIndex(nE, ii, jj, kk);
		if(_dim < 2 && jj != 0)
			continue;
		if(_dim < 3 && kk != 0)
			continue;
		// element state
		// for element 1: composition + state
		Kernel::M_TYPE mt = pME->getMaterial();
		vector<double> comp_el(_materials.size(), 0.);
		const IO::MaterialProperties &element = Domains::global()->PM()->getMaterial(mt);
		float comp = pME->getBAtoms();
		if(comp != 0)
			comp = comp / pME->getVolume() / element._densityCm3 * 1e21; // atm/(nm)^3
		comp_el[_material2idx[mt]] = 1 - comp;
		if(comp != 0)
			comp_el[_material2idx[mt] + 1] = comp;
		feliks::elementState es(_materials.size());
		double initx=0, inity=0, initz=0;
		     //hash     number
		map<string, unsigned> defects;
		setMechanicalState(pME, initx, inity, initz, defects);
		_elstates[i++].setState(comp_el, _pDomain->_pRM->getT(), initx, inity,  initz, defects);     //  temp, eps1, eps2, eps3
	}

	const_cast<model &>(_theAnalysis->getMesh()).updateCurrentState(dofset::tn1);
    
    static bool isFELIKSInitialized=false;
    if ( !isFELIKSInitialized)
    {
        _theAnalysis->solve();
        isFELIKSInitialized = true;
    }
    else
        _theAnalysis->getTheStepSolver().solveStep(logger::mainlog);
    
    
	const_cast<model &>(_theAnalysis->getMesh()).updateCurrentState(dofset::tn1);
	for(unsigned nE = 0; nE < _pDomain->_pMesh->size(); ++nE)
	{
		Kernel::MeshElement * pME = _pDomain->_pMesh->getElement(nE);
		unsigned ii, jj, kk;
		_pDomain->_pMesh->getIndicesFromIndex(nE, ii, jj, kk);
		unsigned i = 0;
		if(_dim == 3)
			i = ii*(_nZ-1)*(_nY-1) + jj*(_nZ-1) + kk;
		if(_dim == 2)
			i = ii*(_nY-1) + jj;
		if(_dim == 1)
			i = ii;

		pME->strain_xx() = _elstates[i].exx;
		pME->strain_yy() = _elstates[i].eyy;
		pME->strain_zz() = _elstates[i].ezz;
		pME->strain_xy() = _elstates[i].exy;
		pME->strain_xz() = _elstates[i].exz;
		pME->strain_yz() = _elstates[i].eyz;

		pME->stress_xx() = _elstates[i].sxx;
		pME->stress_yy() = _elstates[i].syy;
		pME->stress_zz() = _elstates[i].szz;
		pME->stress_xy() = _elstates[i].sxy;
		pME->stress_xz() = _elstates[i].sxz;
		pME->stress_yz() = _elstates[i].syz;
	}
	_theAnalysis->infoSummary(logger::sumf);
	LOWMSG("done.");
}

//physical implementation of the model.
void FELIKSMechanics::setMechanicalState(const Kernel::MeshElement * pME, double &x, double &y, double &z, map<string, unsigned> &defects) const
{
	//first contribution: amorphous expansion.
	x = y = z = Domains::global()->PM()->getMaterial(pME->getMaterial()).amorphousExpansion_;
	//second contribution: number of vacancies.
	const OKMC::Particle *pPart = pME->getFirstPart();
	std::set<OKMC::Defect *> already;
	while(pPart)
	{
		OKMC::Defect *pDef = pPart->getDefect();
		if(pDef->getElement() == pME && already.find(pDef) == already.end())  //the particle contains a defect to this box, and is not accounted already
		{
			already.insert(pDef);
			Kernel::ID id = pDef->getID();
			string name = Domains::global()->PM()->getIDName(id);
			map<string, unsigned>::iterator it = defects.find(name);
			if(it != defects.end())
				it->second++;
			else
				defects[name] = 1;
		}
		pPart = pPart->getNext();
	}
}

/*  collect here a series of initialization tasks that need to be done no matter
 what type of analysis we perform
 */
static void initialize()
{
    // compute machine epsilon
    char   e='E';
	global_macheps = dlamch_(&e, 1);

    // this is called once, and we use it to seed the random numbers
    srand48((long) time(0));
}


}
#else
namespace Mechanics 
{
        FELIKSMechanics::FELIKSMechanics(Kernel::Domain *p) : MechInterface(p) { ERRORMSG("FEM: FELIKS module not compiled!!. You should reconfigure and recompile!"); }
        unsigned FELIKSMechanics::idx2vert(unsigned i, unsigned j, unsigned k) const { return 0; }
        FELIKSMechanics::~FELIKSMechanics() { }
        void FELIKSMechanics::import() const { }
        void FELIKSMechanics::setMechanicalState(const Kernel::MeshElement * pME, double &x, double &y, double &z, std::map<std::string,unsigned> &) const { }
}

#endif
