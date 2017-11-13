/* Author: ignacio.martin@imdea.org
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

#include "Domain.h"
#include "Mesh.h"
#include "RateManager.h"
#include "MeshParam.h"
#include "NUTCreator.h"
#include "electrostatics/Poisson.h"
#include "domains/MCClient.h"
#include "domains/Global.h"
#include "lkmc/LatticeParam.h"
#include "lkmc/LKMCModel.h"
#include "lkmc/Lattice.h"
#include "lkmc/EpiGasParam.h"
#include "io/FileParameters.h"
#include "io/ParameterManager.h"
#include "okmc/MobileParticleParam.h"
#include "okmc/AlloyParam.h"
#include "okmc/InterfaceParam.h"
#include "okmc/ClusterParam.h"
#include "okmc/ClusterReactionParam.h"
#include "okmc/OKMCModel.h"
#include "mechanics/MechInterface.h"

#ifdef _CHARGE_MODEL_
#include "charge/Fermi.h"
#include "charge/FermiParam.h"
#include "charge/ChargeManager.h"
#endif
#include "okmc/ChargedStates.h"

using namespace Electrostatics;
using namespace Kernel;

Domain::Domain(unsigned num, Tcl_Interp *pTcl, const Domains::MCClient *, const Coordinates &m, const Coordinates &M) :
		_rng_dom(Domains::global()->getFileParameters()->getInt("MC/General/random.seed")),
		_um_mechani("Mechanics/General/update"), _um_poisson("MC/Electrostatic/update"), _domain_number(num)
{
	_pMePar = new MeshParam(Domains::global()->PM(), Domains::global()->getFileParameters());
    _pAlPar = new OKMC::AlloyParam(Domains::global()->PM(), Domains::global()->getFileParameters());
    _pMPPar = new OKMC::MobileParticleParam(Domains::global()->PM(), Domains::global()->getFileParameters());

	Coordinates mm = m, MM = M;
	NUTCreator nut(this, mm, MM);
	_pRM = new RateManager(this);
	_pMesh = new Mesh(this, mm, MM,
		nut.getLines(0), nut.getLines(1), nut.getLines(2),
		Domains::global()->client());
#ifndef _CHARGE_MODEL_
		_pSM = new OKMC::ChargedStates(this);
#else
	_pSM = new Charge::ChargeManager(this);
#endif

    _pClPar = new OKMC::ClusterParam(pTcl, Domains::global()->PM(), Domains::global()->getFileParameters());
	const_cast<OKMC::MobileParticleParam *>(_pMPPar)->init(Domains::global()->PM(), Domains::global()->getFileParameters(), _pClPar, _pAlPar, _pSM);
	_pClRePar = new OKMC::ClusterReactionParam(_pClPar, Domains::global()->getFileParameters());
	_pIFPar = new OKMC::InterfaceParam(Domains::global()->PM(), Domains::global()->getFileParameters(),
			const_cast<OKMC::MobileParticleParam *>(_pMPPar), _pClPar, _pMePar);
	for(M_TYPE mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
	{
		_pLaPar[mt] = LKMC::readLatticeParams(pTcl, Domains::global()->PM(), Domains::global()->getFileParameters(), mt);
		_pLat[mt] = LKMC::readLattice(this, _pLaPar[mt], mt);
		_pEGPar[mt] = new LKMC::EpiGasParam;
	}
	Domains::global()->checkInteractions(Domains::global()->PM(), Domains::global()->getFileParameters());

	_pLKMC = 0;
	_pOKMC = 0;

	if (Domains::global()->getFileParameters()->getBool("MC/Electrostatic/load.poisson"))
		_pPoisson = new Poisson(pTcl, this);
	else
		_pPoisson = NULL;

	_pMech = Mechanics::MechInterface::get(this);
}

void Domain::init(bool bFromStream)
{
	_pOKMC = new OKMC::OKMCModel(this);
	_pLKMC = new LKMC::LKMCModel(bFromStream, this);


#ifdef _CHARGE_MODEL_
	Charge::ChargeManager *pCM = dynamic_cast<Charge::ChargeManager *>(_pSM);
	if(pCM)
		pCM->initCharge();
#endif

}

Domain::~Domain()
{
	delete _pMech;
	delete _pLKMC;
	delete _pOKMC;
	delete _pRM;
	delete _pMesh;
	for(M_TYPE mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
	{
		delete _pLat[mt];
		delete _pLaPar[mt];
		delete _pEGPar[mt];
	}
	delete _pSM;
	delete _pIFPar;
	delete _pClRePar;
	delete _pClPar;
    delete _pAlPar;
	delete _pMPPar;
	delete _pMePar;
	delete _pPoisson;
}
