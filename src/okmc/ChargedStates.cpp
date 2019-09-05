/*
 * ChargedStates.cpp
 *
 *  Created on: Aug 23, 2013
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

#include "ChargedStates.h"
#include "domains/Global.h"
#include "io/ParameterManager.h"
#include "kernel/Material.h"
#include "kernel/MeshElement.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/RateManager.h"
#include "MobileParticleParam.h"
#include "Defect.h"

//ChargedStates is mesh element centered instead of particle centered to avoid checking the diffusion barriers
//at every jump

using Kernel::P_TYPE;
using Kernel::MID_STATE;
using Kernel::UNDEFINED_STATE;
using Kernel::UNDEFINED_TYPE;
using Kernel::MAX_STATES;

namespace OKMC {

ChargedStates::ChargedStates(Kernel::Domain *p) : StateManager(p)
{
}

ChargedStates::~ChargedStates()
{
}

unsigned ChargedStates::changeState(Kernel::SubDomain *pSub, P_TYPE pt, const Kernel::MeshElement *pME) const
{
	Kernel::M_TYPE mt = pME->getMaterial();
	const OKMC::MobileParticleParam *pMPPar = _pDomain->_pMPPar;
	if(pMPPar->_bCharged[mt] == false)
		return 0;
    unsigned NStates = Domains::global()->PM()->getStates(mt, pt);
    if(NStates == 1)
    	return 0;
	double eF = pME->bandGap()/2. + pME->electrostaticPotential(); //fermi to valence band
	double scale = pME->bandGap() / Domains::global()->PM()->getMaterial(mt)._Eg0;  //how to scale band gap values
	double probs[Kernel::MAX_STATES];
	const float *levels = pMPPar->_chargeLevel[mt][pt];
	double kT = _pDomain->_pRM->getkT();
	double total = 0;
	for(unsigned i=0; i<NStates; ++i)
	{
		probs[i] = 0;
		int charge = _pDomain->_pMPPar->_state2charge[mt][pt][i];
		switch(charge)
		{
		case 0:
			probs[i] = 1;
			break;
		case 1:
			probs[i] = exp(-(eF - levels[MID_STATE+1]*scale)/kT);
			break;
		case 2:
			probs[i] = exp(-(2*eF - (levels[MID_STATE+2] - levels[MID_STATE+1])*scale)/kT);
			break;
		case 3:
			probs[i] = exp(-(3*eF - (levels[MID_STATE+3] - levels[MID_STATE+2] - levels[MID_STATE+1])*scale)/kT);
			break;
		case -1:
			probs[i] = exp(-(scale*levels[MID_STATE-1] - eF)/kT);
			break;
		case -2:
			probs[i] = exp(-(scale*(levels[MID_STATE-1] + levels[MID_STATE-2]) - 2*eF)/kT);
			break;
		case -3:
			probs[i] = exp(-(scale*(levels[MID_STATE-1] + levels[MID_STATE-2] + levels[MID_STATE-3]) - 3*eF)/kT);
			break;
		}
		total += probs[i];
	}
	//pick one
	float rn = pSub->_rng.rand() * total;
	unsigned lastValid = 0;
	for(unsigned i=0; i<NStates; ++i)
		if(probs[i])
		{
			lastValid = i;
			if (probs[i] > rn)
				return i;
			rn -= probs[i];
		}
	return lastValid;
}

//from pME1 to pME2
bool ChargedStates::checkStateBarrier(Kernel::SubDomain *pSub, P_TYPE pt, unsigned st, const Kernel::MeshElement *pME1, const Kernel::MeshElement *pME2) const
{
	Kernel::M_TYPE mt = pME1->getMaterial();
	if(_pDomain->_pMPPar->_bCharged[mt] == false)
		return true;
	double Eg0 =  Domains::global()->PM()->getMaterial(mt)._Eg0;
	float scale1 = pME1->bandGap() / Eg0;
	float scale2 = pME2->bandGap() / Eg0;
	float eF1   = pME1->bandGap()/2. + pME1->electrostaticPotential();
	float eF2   = pME2->bandGap()/2. + pME2->electrostaticPotential();

	int charge = _pDomain->_pMPPar->_state2charge[mt][pt][st];
	//WARNINGMSG("Partícula " << Domains::global()->PM()->getParticleName(mt, pt) << " tiene carga " << charge);
	float level =_pDomain->_pMPPar->_chargeLevel[mt][pt][st];
	double energy = (eF2 - eF1 + level*(scale1 - scale2))*charge;
	if(energy < 0)
		return true;
	return pSub->_rng.rand() < exp(-energy/_pDomain->_pRM->getkT());
}

double ChargedStates::bindingShift(P_TYPE pt, Kernel::P_POS pos, unsigned st, const Kernel::MeshElement *pME) const
{
	Kernel::M_TYPE mt   = pME->getMaterial();
	const OKMC::MobileParticleParam *pMPPar = _pDomain->_pMPPar;
	if(pMPPar->_bCharged[mt] == false)
		return 0;
	P_TYPE emit, emit2;
	Domains::global()->PM()->getBreakComponents(mt, pt, pos, emit, emit2);
	double Eg0       = Domains::global()->PM()->getMaterial(mt)._Eg0;

	if(emit == UNDEFINED_TYPE || emit2 == UNDEFINED_TYPE)
		return 0;
	float scale = pME->bandGap() / Eg0;
	int imp_charge      = pMPPar->_state2charge[mt][emit][0];
	int my_charge       = pMPPar->_state2charge[mt][pt][st];
	const float *levDop = pMPPar->_chargeLevel[mt][pt];
	const float *levIV  = pMPPar->_chargeLevel[mt][emit2];
	const int diff = my_charge - imp_charge;
	if(diff == 0)
		return 0;
	double correction = 0;
	if(diff > 0) // +1  Bi0 -> B-  I0 -> I+
	{
		int newSt = imp_charge != -1? imp_charge + 1 : imp_charge;
		correction += levDop[MID_STATE + newSt] - levIV[MID_STATE+1];
	}
	if(diff > 1) // +2
		correction += levDop[MID_STATE + imp_charge + 2] - levIV[MID_STATE+2];
	if(diff > 2) // +3
		correction += levDop[MID_STATE + imp_charge + 3] - levIV[MID_STATE+3];
	if(diff < 0) //-1
	{
		int newSt = imp_charge != 1? imp_charge - 1 : imp_charge;
		correction += levIV[MID_STATE-1] - levDop[MID_STATE + newSt];
	}
	if(diff < -1) // -2
		correction += levIV[MID_STATE-2] - levDop[MID_STATE + imp_charge - 2];
	if(diff < -2) // -3
		correction += levIV[MID_STATE-3] - levDop[MID_STATE + imp_charge - 3];
	return correction*scale;
}

unsigned ChargedStates::interactionState (P_TYPE pt1, unsigned st1, P_TYPE pt2, unsigned st2, const Kernel::MeshElement *pME, P_TYPE pt3) const
{
	Kernel::M_TYPE mt   = pME->getMaterial();
	const OKMC::MobileParticleParam *pMPPar = _pDomain->_pMPPar;
	if(_pDomain->_pMPPar->_bCharged[mt] == false)
		return 0;
	int ch1  = pMPPar->_state2charge[mt][pt1][st1];
	int ch2  = pMPPar->_state2charge[mt][pt2][st2];
	int ch3  = ch1 + ch2;
	if(ch3 > int(MID_STATE) || ch3 < -int(MID_STATE))
		return UNDEFINED_STATE;
	return pMPPar->_charge2state[mt][pt3][MID_STATE + ch3];
}

// A_X^j -> A^i + X^j-i
unsigned ChargedStates::breakUpState(P_TYPE pt, Kernel::P_POS pos, unsigned st, Kernel::M_TYPE mt) const
{
	if(_pDomain->_pMPPar->_bCharged[mt] == false)
		return 0;
	P_TYPE emit, emit2;
	Domains::global()->PM()->getBreakComponents(mt, pt, pos, emit, emit2);
	assert(emit != UNDEFINED_TYPE && emit2 != UNDEFINED_TYPE);
	int neutral_charge  = _pDomain->_pMPPar->_state2charge[mt][emit][0];
	int my_charge       = _pDomain->_pMPPar->_state2charge[mt][pt][st];
	return _pDomain->_pMPPar->_charge2state[mt][emit2][MID_STATE + my_charge - neutral_charge];
}

bool ChargedStates::canInteract(P_TYPE pt1, unsigned st1, P_TYPE pt2, unsigned st2, const Kernel::MeshElement *pME) const
{
	Kernel::M_TYPE mt   = pME->getMaterial();
	if(_pDomain->_pMPPar->_bCharged[mt] == false)
		return true;
	//allow all non repulsive reactions.
	int ch1 = _pDomain->_pMPPar->_state2charge[mt][pt1][st1];
	int ch2 = _pDomain->_pMPPar->_state2charge[mt][pt2][st2];

	return ch1*ch2 <= 0;
}


} /* namespace Mechanics */
