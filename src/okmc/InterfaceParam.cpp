/*
 * InterfaceParam.cpp
 *
 *  Created on: Feb 28, 2011
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

#include "InterfaceParam.h"
#include "io/FileParameters.h"
#include "io/ParameterManager.h"
#include "ClusterParam.h"
#include "MobileParticleParam.h"
#include "kernel/MeshParam.h"

using Kernel::M_TYPE;
using Kernel::P_TYPE;
using Kernel::V_TYPE;
using Kernel::MAX_PARTICLES;
using Kernel::P_POS;
using Kernel::NO_POS;

using std::string;
using std::map;
using std::vector;

namespace OKMC {

//initialize everything to 0.
InterfaceParam::InterfaceParamSet::InterfaceParamSet(const OKMC::ClusterParam *pMCPar, M_TYPE mt)
{
	for(int pt=0; pt < MAX_PARTICLES; ++pt)
	{
		_trapMPProb[pt] = 1.;
		_desorptionMPProbL[pt] = _desorptionMPProbH[pt] = _desorptionThreshold[pt] = 0;
		_barrierMP[pt] = 0;
		_interactWithMP[pt] = false;
		_superSat[pt] = 0;
	}

	for(unsigned edtype=0; edtype < pMCPar->defectSize(mt); ++edtype)
	{
		_interactWithMC.push_back(false);
		_trapMCProb.push_back(1.);
		_desorptionMCProb.push_back(0);
	}
	_lambda = 0.5;
}


InterfaceParam::InterfaceParam(const IO::ParameterManager *pPM, const IO::FileParameters * pPar,
		OKMC::MobileParticleParam * pMPPar, const OKMC::ClusterParam *pMCPar, const Kernel::MeshParam *pMePar)
{
	_nMat = pPM->getNMaterials();
	for(M_TYPE mt0=0; mt0 < _nMat; ++mt0)
		for(M_TYPE mt1=0; mt1 < _nMat; ++mt1)
		{
			string interName = getInterfaceName(pPM, mt0, mt1);
			if(pPM->isGasLike(mt1) && pPM->isGasLike(mt0) && mt0 != mt1)
				WARNINGMSG("The interface " << interName << " between two materials without particles is probably wrong");
			if(pPar->specified(interName) && !pPM->isGasLike(mt1))
				_params[mt0][mt1] = new InterfaceParamSet(pMCPar, mt1);
			else
				_params[mt0][mt1] = 0;
		}
	readReactions(pPM, pPar, pMCPar);
	readParameters(pPM, pPar, pMPPar, pMCPar, pMePar);
}

InterfaceParam::~InterfaceParam()
{
	for(M_TYPE mt0=0; mt0 < _nMat; ++mt0)
		for(M_TYPE mt1=0; mt1 < _nMat; ++mt1)
			delete _params[mt0][mt1];
}

string InterfaceParam::getInterfaceName(const IO::ParameterManager *pPM, M_TYPE mt1, M_TYPE mt2) const
{
	string mat1 = pPM->getMaterialName(mt1);
	string mat2 = pPM->getMaterialName(mt2);

	if(mat1 < mat2)
		return mat1 + "_" + mat2;
	else
		return mat2 + "_" + mat1;
}

string InterfaceParam::getTo(const IO::ParameterManager *pPM, const IO::FileParameters * pPar, const std::string &base, M_TYPE mtTo) const
{
	if (pPar->getString(base+"/left") == pPM->getMaterialName(mtTo))
		return "left";
	else
		return "right";
}

void InterfaceParam::readParameters(
		const IO::ParameterManager *pPM,
		const IO::FileParameters * pPar,
		OKMC::MobileParticleParam * pMPPar,
		const OKMC::ClusterParam *pMCPar,
		const Kernel::MeshParam *pMePar)
{
	//It corrects prefactors so that [X] is the concentration written by the user in the input file.
	for(M_TYPE mt0=0; mt0 < pPM->getNMaterials(); ++mt0)
		for(M_TYPE mt1=0; mt1 < pPM->getNMaterials(); ++mt1)
		{
			if(_params[mt0][mt1] == 0)
				continue;
			string interName = getInterfaceName(pPM, mt0, mt1);
			P_TYPE pt1 = Domains::global()->PM()->getMaterial(mt1)._pt[0];
			P_TYPE pt2 = Domains::global()->PM()->getMaterial(mt1)._pt[1];
			double const lambdaNM = pMePar->_lambda[mt1];
			// D = 1/4 * lambda^2 /tau. 1e-7 = cm to nm  4 in 2D
			_params[mt0][mt1]->_lambda = pPar->getFloat(interName + "/Models/lambda");
			double const convFactorMig = 4./(_params[mt0][mt1]->_lambda*_params[mt0][mt1]->_lambda)*1e14;
			for(unsigned fam=0; fam < pPM->getNFamilies(); ++fam)
			{
				if(!pPM->isParticleDefined(fam, mt1))
					continue;
				string baseI = string(interName + "/" + pPM->getFamilyName(fam)); //for interfaces
				string baseB = string(pPM->getMaterialName(mt1) + "/" + pPM->getFamilyName(fam)); //for bulk
				string To    = getTo(pPM, pPar, baseI, mt1);
				// Read all the map<string,float> parameters first
				// left / right parameters
				{
					static string paramNames[3] = { "recombination.length.", "barrier.", "supersaturation." };
					for(unsigned param = 0; param < 3; ++param)
					{
						const bool selfDefect = fam == V_TYPE || fam == pt1 || fam == pt2;
						if(!selfDefect && param == 2) //only for supersaturation
							continue;
						map<string, float> theMap = pPar->getFloatMap(baseI + "/" + paramNames[param] + To);
						for(map<string, float>::iterator it=theMap.begin(); it!=theMap.end(); ++it)
						{
							P_TYPE pt = pPM->getParticleNumber(mt1, it->first);
							unsigned number = pMCPar->getDefectNumber(mt1, it->first);
							if(pt != Kernel::UNDEFINED_TYPE)
							{
								if(pPM->getPPos(pt) == NO_POS)
									ERRORMSG(baseI + "/" + paramNames[param] + To << " " << it->first << " has no position!");
								if(param == 0)  //recombination.length.(Lr)  P(trap) = lambda/Lr.
								{
									if(it->second > lambdaNM)
										_params[mt0][mt1]->_trapMPProb[pt] = lambdaNM/it->second;
									else
									{
										LOWMSG(baseI + "/" + paramNames[param] + To + it->first <<
												" can't be lower than lambda. Reassigning interface trapping probability to 1" );
										_params[mt0][mt1]->_trapMPProb[pt] = 1.;
									}
								}
								else if(param == 1)//barrier
									_params[mt0][mt1]->_barrierMP[pt] = it->second;
								else //supersaturation
									_params[mt0][mt1]->_superSat[pt] = it->second;

								if (!pPM->isParticleDefined(pt, mt1))
									ERRORMSG(baseI << " Particle " << it->first << " not defined for TO material.");

							}
							else if(param == 0 && number != pMCPar->defectSize(mt1))
								_params[mt0][mt1]->_trapMCProb[number] = lambdaNM/it->second;
							else
								ERRORMSG(baseI << " Species " << it->first << " not supported.");
						}
					}
				}
				// not left/right parameters
				{
					static string paramNames[3] = { "desorption.low", "desorption.high", "desorption.threshold" };
					for(unsigned param = 0; param < 3; ++param)
					{
						map<string, float> theMap = pPar->getFloatMap(baseI + "/" + paramNames[param]);
						for(map<string, float>::iterator it=theMap.begin(); it!=theMap.end(); ++it)
						{
							P_TYPE pt = pPM->getParticleNumber(mt1, it->first);
							unsigned number = pMCPar->getDefectNumber(mt1, it->first);
							if(pt != Kernel::UNDEFINED_TYPE)
							{
								if (!pPM->isParticleDefined(pt, mt1) && !pPM->isParticleDefined(pt, mt0))
									ERRORMSG(baseI << "/" << paramNames[param] << " Particle " << it->first << " not defined in any side");
								if(pt == pt1 || pt == pt2 || pt == V_TYPE)
									ERRORMSG("Desorption probability cannot be defined in: "<<baseI + "/" + paramNames[param] << " for " << it->first);
								if(pPM->getPPos(pt) == Kernel::NO_POS && (param == 0 || param == 1))
									ERRORMSG(baseI << ". Parameter " << paramNames[param] << " not defined for " << it ->first << ". Particle position needed.");
								P_TYPE family = Domains::global()->PM()->getFamily(pt);
								bool selfDefect = family == V_TYPE || family == pt1 || family == pt2;
								switch(param)
								{
								case 0:
									_params[mt0][mt1]->_desorptionMPProbL[pt] = it->second;
									break;
								case 1:
									_params[mt0][mt1]->_desorptionMPProbH[pt] = it->second;
									break;
								case 2:
									pt = pPM->getFamily(pt);
									_params[mt0][mt1]->_desorptionThreshold[pt] = it->second;
									if (selfDefect)
										ERRORMSG(baseI << ". Parameter " << paramNames[param] << " not defined for " << it ->first);
									break;
								default:
									ERRORMSG("Case not implemented in InterfaceParam!");
								}
							}
							else if(number != pMCPar->defectSize(mt1))
							{
								switch(param)
								{
								case 1:
									_params[mt0][mt1]->_desorptionMCProb[number] = it->second;
									break;
								default:
									ERRORMSG("Parameter " << paramNames[param] << " not defined for " << it->first);
								}
							}
							else
								WARNINGMSG(baseI  + "/" + paramNames[param] << " Species " << it->first << " not supported in " << pPM->getMaterialName(mt1));
						}
					}
				}
				//Read formation energy at the surface
				if(fam != pt1 && fam != pt2 && fam != V_TYPE)
					if(pPM->isParticleDefined(fam, mt0) || pPM->isParticleDefined(fam, mt1)) //pure impurities
					{
						_params[mt0][mt1]->_formation[fam] = pPar->getArrheniusAlloys(baseI + "/" + pPM->getParticleName(mt1, fam) + "(formation)")._v[0];
						_params[mt0][mt1]->_migration[fam] = pPar->getArrheniusAlloys(baseI + "/" + pPM->getParticleName(mt1, fam) + "(migration)")._v[0]*convFactorMig;
					}
			}
			//assign rates
			for(P_TYPE pt=pPM->getNFamilies(); pt < pPM->getNParticles(); ++pt)
			{
				if (!pPM->isParticleDefined(pt, mt1))
					continue;
				P_TYPE fam = Domains::global()->PM()->getFamily(pt);
				bool selfDefect = fam == V_TYPE || fam == pt1 || fam == pt2;
				if (!selfDefect)
				{
					//nu_emiss = [A_vol] nu_m lambda / (6 [A_surf]
					//Prefactor: nu_m0 * lambda * CA/CARef / (6 * CA_surf/CARef)
					//energy:  (Ef - EfRef) + Em - (Ef_surf - Ef_Ref)
					// have lambda / (6* CA_surf/CARef), ener = -Ef_surf
					_params[mt0][mt1]->_arrEmitMP[pt]  = IO::ArrheniusAlloys(lambdaNM*1e-7 / _params[mt0][mt1]->_formation[fam]._pref / 6.,
							-_params[mt0][mt1]->_formation[fam]._ener);
					_params[mt0][mt1]->_arrEmitMP[pt] *= pMPPar->_orig[mt1][pt][0][0];
					_params[mt0][mt1]->_arrEmitMP[pt].addEner(_params[mt0][mt1]->_barrierMP[pt]); //add barrier
					//have nu_m0 * lambda / (6* CA_surf/CARef), ener = Ef + Em -Ef_surf
					_params[mt0][mt1]->_arrEmitMP[pt] *= pMPPar->_form[mt1][pt];
				}
				else
				{
					//emissions
					_params[mt0][mt1]->_arrEmitMP[pt] = pMPPar->_form[mt1][pt];
					_params[mt0][mt1]->_arrEmitMP[pt] *= pMPPar->_orig[mt1][pt][0][0] / (lambdaNM*1e-7);
					_params[mt0][mt1]->_arrEmitMP[pt].addEner(_params[mt0][mt1]->_barrierMP[pt]);
				}
				if(_params[mt0][mt1]->_interactWithMP[pt] == false)
					_params[mt0][mt1]->_arrEmitMP[pt] = IO::ArrheniusAlloys();
			}
		}
}

void InterfaceParam::readReactions(const IO::ParameterManager *pPM, const IO::FileParameters * pPar,
		const OKMC::ClusterParam *pMCPar)
{
	for(M_TYPE mt0=0; mt0 < pPM->getNMaterials(); ++mt0)
		for(M_TYPE mt1=0; mt1 < pPM->getNMaterials(); ++mt1)
		{
			if(_params[mt0][mt1] == 0)
				continue;

			string mat0 = pPM->getMaterialName(mt0);
			string mat1 = pPM->getMaterialName(mt1);

			string base = mat1 + "/Models/interactions";
			IO::array<string, string> theMap = pPar->getArray(base);
			for(IO::array<string, string>::iterator it=theMap.begin(); it!=theMap.end(); ++it)
			{
				vector<string> tokens;
				IO::ParameterManager::getTokens(it->first, '+', tokens);
				if(tokens.size() != 2)
					ERRORMSG(base << " An interaction needs only two tokens: " << it->first);
				if(tokens[1] != mat0)  //not this interface
					continue;
				if(it->second == "false")
					continue;
				unsigned edtype = pMCPar->getDefectNumber(mt1, tokens[0]);
				P_TYPE pt = pPM->getParticleNumber(mt1, tokens[0]);
				if(edtype < pMCPar->defectSize(mt1))
				{
					_params[mt0][mt1]->_interactWithMC[edtype] = true;
					Domains::global()->setInteraction(mt1, it->first);
				}
				else if(pt != Kernel::UNDEFINED_TYPE && pPM->getPPos(pt) != NO_POS)
				{
					if(!pPM->isParticleDefined(pt, mt1))
						ERRORMSG(base << " Particle " << tokens[0] << " not defined in material");
					_params[mt0][mt1]->_interactWithMP[pt] = true;
					Domains::global()->setInteraction(mt1, it->first);
				}
				else
					WARNINGMSG(base << " Could not understand reaction " << it->first);

				MEDMSG("Interaction: " << tokens[0] << " + " << "interface in " << pPM->getMaterialName(mt1));
			}
		}
}

}  //end namespace


