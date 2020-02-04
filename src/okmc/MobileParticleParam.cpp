/*
 * MobileParticleParam.cpp
 *
 *  Created on: Feb 24, 2011
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

#include "MobileParticleParam.h"
#include "ClusterParam.h"
#include "Defect.h"
#include "io/FileParameters.h"
#include "io/ParameterManager.h"
#include "kernel/Domain.h"
#include "kernel/MeshElement.h"
#include "kernel/StateManager.h"
#include <cstdlib>

using IO::Arrhenius;
using IO::ArrheniusAlloys;
using Kernel::M_TYPE;
using Kernel::P_TYPE;
using Kernel::V_TYPE;
using Kernel::MID_STATE;
using Kernel::UNDEFINED_STATE;
using Kernel::P_POS;
using Kernel::POS_0;
using Kernel::POS_1;
using Kernel::POS_I;
using Kernel::POS_V;
using Kernel::NO_POS;
using Kernel::MAX_PARTICLES;
using Kernel::UNDEFINED_TYPE;

using std::map;
using std::string;
using std::vector;
using std::stringstream;

namespace OKMC {

MobileParticleParam::MobileParticleParam(const IO::ParameterManager *pPM, const IO::FileParameters * pPar)
{
	//Charge properties
	for(M_TYPE mt = 0; mt < pPM->getNMaterials(); ++mt)
	{
		_bCharged[mt] = false;
		int chargeModel = -1; //uninitialized
		std::stringstream tempMsg;
		for(P_TYPE pt = 0; pt < pPM->getNParticles(); ++pt)
		{
			_stateEnergy[mt][pt]     = 0;
			_stateDegeneracy[mt][pt] = 0;
			_orbitalRadius[mt][pt]   = 0;
			_mapToGrid[mt][pt]       = false;

			if(!pPM->isParticleDefined(pt, mt))
				continue;
			stringstream base;
			base << pPM->getMaterialName(mt) << "/";
			if(pt >= pPM->getNFamilies())
			{
				Kernel::P_TYPE fam = pPM->getFamily(pt);
				base << pPM->getFamilyName(fam) << "/" << pPM->getParticleName(mt, pt);
				_bCharged[mt] = pPar->specified(base.str()+ "(state.charge)");
				tempMsg << base.str()+"(state.charge)" << " is " << (_bCharged[mt]? "" : "not ") << "defined\n";
				if(chargeModel == -1) //checked that chargeModel is consistent in the whole material.
					chargeModel = _bCharged[mt]? 1:0;
				else if((chargeModel == 1 && !_bCharged[mt]) || (chargeModel == 0 && _bCharged[mt]))
					ERRORMSG(base.str() << ". Inconsistency in 'state.charge' definition. Define ALL or NONE.\n"
						<< tempMsg.str());
			}
			for(unsigned st=0; st<Kernel::MAX_STATES; ++st)
				_charge2state[mt][pt][st] = UNDEFINED_STATE;
			if(_bCharged[mt])
			{
				P_TYPE fam = pPM->getFamily(pt);
				P_POS  pos = pPM->getPPos(pt);
				if(fam != V_TYPE && fam != pPM->getMaterial(mt)._pt[0] && fam != pPM->getMaterial(mt)._pt[1] &&
					 pos != POS_I && pos != POS_V)
				{
					_stateEnergy[mt][pt]     = pPar->getFloat(base.str() + "(state.energy)");
					_stateDegeneracy[mt][pt] = pPar->getFloat(base.str() + "(state.degeneracy)");
					if (pPar->specified(base.str() + "(map.to.grid)"))
					{
						_mapToGrid[mt][pt]       = pPar->getBool (base.str() + "(map.to.grid)");
						_orbitalRadius[mt][pt]   = pPar->getFloat(base.str() + "(orbital.radius)");
					}
				}
				if(pt >= pPM->getNFamilies()) //only for "real" particles with position
				{
					map<string, int> chargeMap = pPar->getIntMap(base.str() + "(state.charge)");
					int max_charge = -20;
					for(unsigned st=0; st<pPM->getStates(mt, pt); ++st)
					{
						string stName = pPM->getParticleStateName(mt, pt, st);
						if(chargeMap.find(stName) == chargeMap.end())
							ERRORMSG(base.str() << ". State '" << stName << "' not found in map");
						int charge = chargeMap[stName];
						_state2charge[mt][pt][st] = charge;
						if(st == 0 && charge != 0 && pPM->getStates(mt, pt) > 1)
							ERRORMSG("When defining multiple charges, the first one must be neutral.");
						_charge2state[mt][pt][MID_STATE + charge] = st;
						if(charge > max_charge)
							max_charge = charge;
					}
					if(max_charge > 3)
						ERRORMSG(base.str() << ". Charge state bigger than 3 defined!");
					for(int ch=-int(MID_STATE); ch < max_charge; ++ch)
					{
						unsigned st = _charge2state[mt][pt][MID_STATE+ch];
						if(!pPM->isPureImpurity(mt, pt) && st != UNDEFINED_STATE) //everything but pure impurities
						{
							stringstream level;
							level << base.str() << "(e(" << ch << "," << ch+1 << "))";
							_chargeLevel[mt][pt][MID_STATE + (ch >= 0? ch+1:ch)] = pPar->getFloat(level.str());
						}
					}
				}
			}
		}
	}
}

void MobileParticleParam::init(const IO::ParameterManager *pPM, const IO::FileParameters * pPar,
		const OKMC::ClusterParam *pMC, const OKMC::AlloyParam *pAP, const Kernel::StateManager *pSM)
{
	//init
	readReactions(pPM, pPar, pMC);

	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)  //read 1D or 3D diffusion.
		for(P_TYPE pt=pPM->getNFamilies(); pt < pPM->getNParticles(); ++pt)
		{
			if(!pPM->isParticleDefined(pt, mt))
				continue;
			for(unsigned st=0; st<pPM->getStates(mt, pt); ++st)
			{
				stringstream ss;
				ss << pPM->getMaterialName(mt) << "/";
				Kernel::P_TYPE fam = pPM->getFamily(pt);
				ss << pPM->getFamilyName(fam) << "/" << pPM->getParticleName(mt, pt);

				if(pPM->getStateName(mt, pt, st).size())
					ss << "_" << pPM->getStateName(mt, pt, st);
				_oneD[mt][pt][st] = pPar->specified(ss.str() + "(axis)");
				if(_oneD[mt][pt][st])
				{
					_axes[mt][pt][st] = pPar->getCoordinates(ss.str() + "(axis)");
					if (_axes[mt][pt][st].abs() == 0) //(0,0,0) == 3D
						_oneD[mt][pt][st] = false;
				}
				else
					_axes[mt][pt][st] = Kernel::Coordinates(0,0,0);
			}
		}
	// Mig Br1 Br2 FT_I FT_V Update || Br1 Br2 FT_I FT_V
	// 0    1   2    3    4    5         6  7    8    9
	//diffusion and update
	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
	{
		float lambda = pPar->getFloat(pPM->getMaterialName(mt) + "/Models/lambda");
		// D = 1/(2*dim) * lambda² /tau. 1e-7 = cm to nm
		//convFactor = 2*dim/(lambda*lambda)*1e14;
		//1D migration. Reads the I(axis) and sets the convFactor properly later.
		for(P_TYPE pt=pPM->getNFamilies(); pt < pPM->getNParticles(); ++pt)
		{
			if(!pPM->isParticleDefined(pt, mt))
				continue;
			stringstream base;
			base << pPM->getMaterialName(mt) << "/";
			Kernel::P_TYPE fam = pPM->getFamily(pt);
			base << pPM->getFamilyName(fam) << "/" << pPM->getParticleName(mt, pt);

			if(pPM->getStates(mt, pt) > 1)  //update
			{
				double convFactor = 6./(lambda*lambda)*1e14;
				_orig[mt][pt][0][5] = pPar->getArrheniusAlloys(base.str() + "(update)");
#ifdef NUMODEL
				_arr [mt][pt][0][10] = _orig[mt][pt][0][5]* convFactor;
#else
				_arr [mt][pt][0][5] = _orig[mt][pt][0][5]* convFactor;
#endif
			}
			for(unsigned st=0; st<pPM->getStates(mt, pt); ++st)  //migration
			{
				if(pPM->getStates(mt, pt) > Kernel::MAX_STATES)
				{
					for(unsigned stCheck=0; stCheck<pPM->getStates(mt, pt); ++stCheck)
						LOWMSG(pPM->getStateName(mt, pt, stCheck) << std::endl);
					ERRORMSG("more states defined than allowed for: " <<pPM->getParticleName(mt, pt)
							<< " in material '" << pPM->getMaterialName(mt) <<std::endl);
				}
				if(!pPM->isStateDefined(mt, pt, st))
					continue;
				stringstream ss;
				ss << base.str();
				if(pPM->getStateName(mt, pt, st).size())
					ss << "_" << pPM->getStateName(mt, pt, st);
				double convFactorMig = _oneD[mt][pt][st]? 2./(lambda*lambda)*1e14 : 6./(lambda*lambda)*1e14;
				_orig[mt][pt][st][0] = pPar->getArrheniusAlloys(ss.str() + "(migration)");
#ifdef NUMODEL
				for (unsigned i = 0; i < 6; i++)
				{
					_arr[mt][pt][st][i] = _orig[mt][pt][st][0] * convFactorMig;
					Domains::global()->PM()->setMigrationRate(mt, st, pt, _arr[mt][pt][st][i]);
				}
#else
				_arr[mt][pt][st][0] = _orig[mt][pt][st][0] * convFactorMig;
				Domains::global()->PM()->setMigrationRate(mt, st, pt, _arr[mt][pt][st][0]);
#endif
			}
		}
	}
	//formations
	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
		for(P_TYPE part=pPM->getNFamilies(); part < pPM->getNParticles(); ++part)
		{
			if(!pPM->isParticleDefined(part, mt))
				continue;
			stringstream ss;
			P_TYPE fam = Domains::global()->PM()->getFamily(part);
			ss << "" << pPM->getMaterialName(mt) << "/" << pPM->getFamilyName(fam) << "/" << pPM->getParticleName(mt, part);
			_form[mt][part] = pPar->getArrheniusAlloys(ss.str() + "(formation)");
		}
	//Calculation of break-up energies.
	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
	{
		for(P_TYPE part=pPM->getNFamilies() ; part < pPM->getNParticles(); ++part)
		{
			if(pPM->isParticleDefined(part, mt) == false)
				continue;
			for(unsigned pos = POS_0; pos <= POS_1; ++pos)
			{
				P_TYPE emit, emit2;
				pPM->getBreakComponents(mt, part, P_POS(pos), emit, emit2);
				if(emit == UNDEFINED_TYPE || emit2 == UNDEFINED_TYPE)
					continue;
				if(!_canInteract[mt][emit][emit2])
					continue;
				const float l = _interact_radius[mt][emit][emit2];
				const double captVol = 4./3.*M_PI * l*l*l *1e-21;
				for(unsigned st=0; st<pPM->getStates(mt, part); ++st)
				{
					_breakUp[mt][part][st]._emit[0] = emit;
					_breakUp[mt][part][st]._emit[1] = emit2;
					if(_interact[mt][emit][emit2] == Kernel::Event::CLUSTER)
						_breakUp[mt][part][st]._isCluster = true;
					unsigned delta = (pos == POS_0? 0:1);
					unsigned ivSt = pSM->breakUpState(part, Kernel::P_POS(pos), st, mt);
					if(ivSt != UNDEFINED_STATE && emit != part && emit != UNDEFINED_TYPE)
					{
						if(!_form[mt][part].canDivide())
							ERRORMSG(pPM->getMaterialName(mt) << "/" << pPM->getParticleName(mt, part) << " formation prefactor cannot be zero");
						_arr[mt][part][st][1+delta] = _form[mt][emit] / _form[mt][part] * _form[mt][emit2] * _arr[mt][emit2][ivSt][0] * captVol;
						_arr[mt][part][st][6+delta] = _form[mt][emit] / _form[mt][part] * _form[mt][emit2] * _arr[mt][emit][st][0] * captVol;
					}
				}
			}
		}
	}
	//F-T emission for particles
	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
	{
		const bool binary = pPM->getMaterial(mt)._binary;
		for(P_TYPE part=pPM->getNFamilies(); part < pPM->getNParticles(); ++part)
			for(char iorv=0; iorv < 2; ++iorv)  // dopant dop emiting iorv
			{
				if(pPM->isParticleDefined(part, mt) == false)
					continue;
				P_POS pos = pPM->getPPos(part);
				P_POS otherPos = (pos == POS_0 && binary)? POS_1 : POS_0;
				P_TYPE fam = pPM->getFamily(part);
				if(pos == POS_I || pos == POS_V || fam == V_TYPE) //only substitutionals
					continue;
				P_TYPE emit1 = pPM->getParticle(mt, fam,  (iorv == Kernel::IV_I? POS_V:POS_I));
				P_TYPE emit2 = pPM->iorv2pt    (mt, iorv, (iorv == Kernel::IV_I? otherPos: pos));
				if(emit1 == UNDEFINED_TYPE || emit2 == UNDEFINED_TYPE)
					continue;
				//state 0?
				const float l = _interact_radius[mt][emit1][emit2];
				const double captVol = 4./3.*M_PI * l*l*l * 1e-21;
				if(!_canInteract[mt][emit1][emit2]) //mic. reversib.
					_arr[mt][part][0][3+iorv] = _arr[mt][part][0][8+iorv] = ArrheniusAlloys();
				else
				{
					if(!_form[mt][part].canDivide())
						ERRORMSG(pPM->getMaterialName(mt) << "/" << pPM->getParticleName(mt, part) << " formation prefactor cannot be zero");
					_arr[mt][part][0][3+iorv]= _form[mt][emit1]/_form[mt][part]*_form[mt][emit2]*_arr[mt][emit1][0][0]*captVol;
					_arr[mt][part][0][8+iorv]= _form[mt][emit1]/_form[mt][part]*_form[mt][emit2]*_arr[mt][emit2][0][0]*captVol;
				}
			}
	}	
}

void MobileParticleParam::readReactions(const IO::ParameterManager *pPM, const IO::FileParameters * pPar,
		const OKMC::ClusterParam *pMC)
{
	//init
	for(P_TYPE pt1=0; pt1 < MAX_PARTICLES; ++pt1)
		for(P_TYPE pt2=0; pt2 < MAX_PARTICLES; ++pt2)
		{
			_int_result[pt1][pt2] = UNDEFINED_TYPE;
			for(M_TYPE mt=0; mt < Kernel::MAX_MATERIALS; ++mt)
			{
				_interact[mt][pt1][pt2] = Kernel::Event::EMPTY;
				_interact_radius[mt][pt1][pt2] = _interact_radiusSq[mt][pt1][pt2] = 0;
				_canInteract[mt][pt1][pt2] = false;
			}
		}
	//read
	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
	{
		const float lambda = pPar->getFloat(pPM->getMaterialName(mt) + "/Models/lambda");
		stringstream sss;
		sss << pPM->getMaterialName(mt) << "/Models/interactions";
		IO::array<string, string> theMap = pPar->getArray(sss.str().c_str());
		for(IO::array<string, string>::iterator it=theMap.begin(); it!=theMap.end(); ++it)
		{
			if(it->second == "false")
				continue;
			vector<string> tokens, tokens2;
			IO::ParameterManager::getTokens(it->first, '+', tokens);
			IO::ParameterManager::getTokens(it->second, ',', tokens2);
			if(tokens.size() != 2)
				ERRORMSG(sss.str() << " An interaction needs only two tokens: " << it->first);
			P_TYPE reactants[2] = { pPM->getParticleNumber(mt, tokens[0]), pPM->getParticleNumber(mt, tokens[1]) };
			if(reactants[0] == UNDEFINED_TYPE || reactants[1] == UNDEFINED_TYPE) //it could be an Interface, or an extended defect...
				continue;
			if(pPM->getPPos(reactants[0]) == NO_POS || pPM->getPPos(reactants[1]) == NO_POS)
				ERRORMSG(sss.str() << " Some particle position " << it->first << " is not defined here.");
			if(!pPM->isParticleDefined(reactants[0], mt) || !pPM->isParticleDefined(reactants[1], mt))
				ERRORMSG(sss.str() << " Some particle " << it->first << " is not defined here.");
			if(pPM->getPPos(reactants[0]) == NO_POS || pPM->getPPos(reactants[1]) == NO_POS)
				ERRORMSG(sss.str() << " Some particle " << it->first << " do not have a valid position.");
			P_TYPE results[MAX_PARTICLES][MAX_PARTICLES];
			bool   iv     [MAX_PARTICLES][MAX_PARTICLES];
			fillResults(mt, results, iv);
			if(iv[reactants[0]][reactants[1]] && tokens2[0] == "0")
			{
				if(tokens2.size() != 1 && tokens2.size() != 2)
					ERRORMSG(sss.str() << " Interaction result syntax error. " << it->second);
				float capRad = tokens2.size() == 2 ? atof(tokens2[1].c_str()) : lambda;
				float capRadSq = capRad * capRad;
				Domains::global()->setInteraction(mt, it->first);
				_interact[mt][reactants[0]][reactants[1]] = Kernel::Event::EMPTY;
				_interact[mt][reactants[1]][reactants[0]] = Kernel::Event::EMPTY;
				_canInteract[mt][reactants[0]][reactants[1]] = true;
				_canInteract[mt][reactants[1]][reactants[0]] = true;
				_interact_radius[mt][reactants[0]][reactants[1]] = capRad;
				_interact_radius[mt][reactants[1]][reactants[0]] = capRad;
				Domains::global()->PM()->setMaximumCaptureRadius(mt, capRad);
				_interact_radiusSq[mt][reactants[0]][reactants[1]] = capRadSq;
				_interact_radiusSq[mt][reactants[1]][reactants[0]] = capRadSq;
				MEDMSG(sss.str().c_str() << " " << tokens[0] << " + " << tokens[1]);
				continue;
			}

			//MobileParticle
			if(results[reactants[0]][reactants[1]] != UNDEFINED_TYPE &&
					(tokens2[0] == "true" || tokens2[0] == "sink"))
			{
				if(tokens2.size() != 1 && tokens2.size() != 2)
					ERRORMSG(sss.str() << " Interaction result syntax error. " << it->second);
				float capRad = tokens2.size() == 2 ? atof(tokens2[1].c_str()) : lambda;
				float capRadSq = capRad * capRad;
				_canInteract[mt][reactants[0]][reactants[1]] = true;
				_canInteract[mt][reactants[1]][reactants[0]] = true;
				if(tokens2[0] == "true")
				{
					_interact[mt][reactants[0]][reactants[1]] = Kernel::Event::MOBILEPARTICLE;
					_interact[mt][reactants[1]][reactants[0]] = Kernel::Event::MOBILEPARTICLE;
				}
				else
				{
					_interact[mt][reactants[0]][reactants[1]] = Kernel::Event::SINK;
					_interact[mt][reactants[1]][reactants[0]] = Kernel::Event::SINK;
				}
				_interact_radius[mt][reactants[0]][reactants[1]] = capRad;
				_interact_radius[mt][reactants[1]][reactants[0]] = capRad;
				Domains::global()->PM()->setMaximumCaptureRadius(mt, capRad);
				_interact_radiusSq[mt][reactants[0]][reactants[1]] = capRadSq;
				_interact_radiusSq[mt][reactants[1]][reactants[0]] = capRadSq;
				Domains::global()->setInteraction(mt, it->first);
				_int_result[reactants[0]][reactants[1]] = results[reactants[0]][reactants[1]];
				_int_result[reactants[1]][reactants[0]] = results[reactants[0]][reactants[1]];
				MEDMSG(sss.str().c_str() << " " << tokens[0] << " + " << tokens[1]);
				continue;
			}

			//Clusters is everything else
			if(tokens2.size() != 2 && tokens2.size() != 3)
				ERRORMSG("Cluster interaction needs the probability to form the cluster (like ICluster,1) '"
						<< sss.str() << " " << it->first << "' -> '" << it->second << "'");

			unsigned edType = pMC->getDefectNumber(mt, tokens2[0]);
			float capRad = tokens2.size() == 3 ? atof(tokens2[2].c_str()) : lambda;
			float capRadSq = capRad * capRad;
			if(edType >= pMC->defectSize(mt))
				ERRORMSG(pPM->getMaterialName(mt) + "/Models/interactions: Extended defect " << tokens2[0] << " not recognized.");
			_mctype[mt][reactants[0]][reactants[1]][edType] = std::atof(tokens2[1].c_str());
			_mctype[mt][reactants[1]][reactants[0]][edType] = std::atof(tokens2[1].c_str());
			_interact[mt][reactants[0]][reactants[1]] = Kernel::Event::CLUSTER;
			_interact[mt][reactants[1]][reactants[0]] = Kernel::Event::CLUSTER;
			_interact_radius[mt][reactants[0]][reactants[1]] = capRad;
			_interact_radius[mt][reactants[1]][reactants[0]] = capRad;
			Domains::global()->PM()->setMaximumCaptureRadius(mt, capRad);
			_interact_radiusSq[mt][reactants[0]][reactants[1]] = capRadSq;
			_interact_radiusSq[mt][reactants[1]][reactants[0]] = capRadSq;
			Domains::global()->setInteraction(mt, it->first);

			_canInteract[mt][reactants[0]][reactants[1]] = true;
			_canInteract[mt][reactants[1]][reactants[0]] = true;

			MEDMSG(sss.str().c_str() << " " << tokens[0] << " + " << tokens[1]);
		}
	}
}

void MobileParticleParam::fillResults(M_TYPE mt,
		P_TYPE pds[MAX_PARTICLES][MAX_PARTICLES],
		bool iv[MAX_PARTICLES][MAX_PARTICLES]) const
{
	for(unsigned i=0; i < MAX_PARTICLES; ++i)
		for(unsigned j=0; j< MAX_PARTICLES; ++j)
		{
			pds[i][j] = UNDEFINED_TYPE;
			iv[i][j] = false;
		}
	if(Domains::global()->PM()->getMaterial(mt)._binary)
	{
		const P_TYPE Sii = Domains::global()->PM()->iorv2pt(mt, Kernel::IV_I, POS_0);
		const P_TYPE Ci  = Domains::global()->PM()->iorv2pt(mt, Kernel::IV_I, POS_1);
		const P_TYPE VSi = Domains::global()->PM()->iorv2pt(mt, Kernel::IV_V, POS_0);
		const P_TYPE VC  = Domains::global()->PM()->iorv2pt(mt, Kernel::IV_V, POS_1);
		const P_TYPE SiC = Domains::global()->PM()->getAntisite(mt, POS_0);
		const P_TYPE CSi = Domains::global()->PM()->getAntisite(mt, POS_1);
		const P_TYPE mt1 = Domains::global()->PM()->getMaterial(mt)._pt[0];
		const P_TYPE mt2 = Domains::global()->PM()->getMaterial(mt)._pt[1];

		iv[Sii][VSi] = iv[VSi][Sii] = true;
		iv[SiC][CSi] = iv[CSi][SiC] = true;
		iv[Ci][VC]   = iv[VC][Ci]   = true;

		pds[Sii][CSi] = Ci;
		pds[Sii][VC]  = SiC;
		pds[SiC][Ci]  = Sii;
		pds[SiC][VSi] = VC;
		pds[Ci] [SiC] = Sii;
		pds[Ci] [VSi] = CSi;
		pds[CSi][Sii] = Ci;
		pds[CSi][VC]  = VSi;
		pds[VSi][Ci]  = CSi;
		pds[VSi][SiC] = VC;
		pds[VC][Sii]  = SiC;
		pds[VC][CSi]  = VSi;

		for(unsigned family=1; family < Domains::global()->PM()->getNFamilies(); ++family)
		{
			if(family == mt1 || family == mt2)
				continue;
			//pure "impurities" now
			P_TYPE Ali  = Domains::global()->PM()->getParticle(mt, family, POS_I);
			P_TYPE AlC  = Domains::global()->PM()->getParticle(mt, family, POS_1);
			P_TYPE AlSi = Domains::global()->PM()->getParticle(mt, family, POS_0);
			P_TYPE AlV  = Domains::global()->PM()->getParticle(mt, family, POS_V);
			pds[Ali][VSi] = pds[VSi][Ali] = AlSi;
			pds[Ali][VC]  = pds[VC][Ali]  = AlC;
			pds[AlC][Ci]  = pds[Ci][AlC]  = Ali;
			pds[AlC][CSi] = pds[CSi][AlC] = AlSi;
			pds[AlC][VSi] = pds[VSi][AlC] = AlV;
			pds[AlSi][Sii]= pds[Sii][AlSi]= Ali;
			pds[AlSi][SiC]= pds[SiC][AlSi]= AlC;
			pds[AlSi][VC] = pds[VC][AlSi] = AlV;
			pds[AlV][Sii] = pds[Sii][AlV] = AlC;
			pds[AlV][Ci]  = pds[Ci][AlV]  = AlSi;
			pds[AlV][CSi] = pds[CSi][AlV] = AlV;
		}
	}
	else
	{
		const P_TYPE Sii = Domains::global()->PM()->iorv2pt(mt, Kernel::IV_I, POS_0);
		const P_TYPE VSi = Domains::global()->PM()->iorv2pt(mt, Kernel::IV_V, POS_0);
		const P_TYPE mt1 = Domains::global()->PM()->getMaterial(mt)._pt[0];

		iv[Sii][VSi] = iv[VSi][Sii] = true;

		for(unsigned family=1; family < Domains::global()->PM()->getNFamilies(); ++family)
		{
			if(family == mt1)
				continue;
			//pure "impurities" now
			P_TYPE Ali  = Domains::global()->PM()->getParticle(mt, family, POS_I);
			P_TYPE AlSi = Domains::global()->PM()->getParticle(mt, family, POS_0);
			P_TYPE AlV  = Domains::global()->PM()->getParticle(mt, family, POS_V);
			pds[Ali][VSi] = pds[VSi][Ali] = AlSi;
			pds[AlSi][Sii]= pds[Sii][AlSi]= Ali;
			pds[AlSi][VSi] = pds[VSi][AlSi] = AlV;
			pds[AlV][Sii] = pds[Sii][AlV] = AlSi;
		}
	}
}

}


