/*
 * ClusterParam.cpp
 *
 *  Created on: May 30, 2011
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

#include "ClusterParam.h"
#include "Cluster.h"
#include "io/ParameterManager.h"
#include "io/FileParameters.h"
#include "EDTypeDisk.h"
#include "EDTypeIrregular.h"
#include "EDTypePlane311.h"
#include "EDTypeVoid.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/MeshElement.h"
#include "okmc/MobileParticleParam.h"
#include "okmc/Defect.h"
#include <cmath>
#include <fnmatch.h>

using std::map;
using Kernel::M_TYPE;
using Kernel::ID;
using Kernel::P_TYPE;
using Kernel::V_TYPE;
using Kernel::UNDEFINED_TYPE;
using Kernel::P_POS;
using Kernel::NO_POS;
using Kernel::POS_0;
using Kernel::POS_1;
using Kernel::POS_V;
using std::vector;
using std::string;

namespace OKMC {

ClusterParam::ClusterParam(Tcl_Interp *pTcl, const IO::ParameterManager *pPM, const IO::FileParameters * pPar) : _pPM(pPM)
{
	readParameters(pTcl, pPM, pPar);
	readInteractions(pPar);
}

void ClusterParam::readParameters(Tcl_Interp *pTcl, const IO::ParameterManager *pPM, const IO::FileParameters * pPar)
{
	LOWMSG2("Starting clusters: ");
	for(M_TYPE mt=0; mt < _pPM->getNMaterials(); ++mt)
	{
		LOWMSG2(" " << _pPM->getMaterialName(mt) + "(");
		map<string,bool> theArray = pPar->getBoolMap(_pPM->getMaterialName(mt) + "/Models/defined");
		const float lambda = pPar->getFloat(_pPM->getMaterialName(mt) + "/Models/lambda");
		// D = 1/6 * lambda^2 /tau. 1e-7 = cm to nm
		double const convFactor = 6./(lambda*lambda)*1e14;
		for(map<string,bool>::iterator it=theArray.begin(); it!=theArray.end(); ++it)
		{
			if(it->second == false)
				continue;
			LOWMSG2(it->first << " ");
			string base_param = pPM->getMaterialName(mt) + '/' + it->first;
			string aspect = pPar->getString(base_param + "/shape");
			EDType * edtype = 0;
			if(aspect == "disk")
			{
				edtype = new EDTypeDisk();
				edtype->_densityNm = 1e-14 * pPar->getFloat(base_param+ "/density.cm2");
				EDTypeDisk *theDisk = static_cast< EDTypeDisk *>(edtype);
				theDisk->_ratio = pPar->getFloat(base_param+"/axes.ratio");
			}
			else if(aspect == "void")
			{
				edtype = new EDTypeVoid();
				edtype->_densityNm = 1e-21 * pPar->getFloat(base_param+ "/density.cm3");
			}
			else if(aspect == "plane311")
			{
				edtype = new EDTypePlane311();
				edtype->_densityNm = 1e-14 * pPar->getFloat(base_param+ "/density.cm2");
			}
			else if (aspect == "irregular")
			{
				edtype = new EDTypeIrregular();
				edtype->_densityNm = 1e-21 * pPar->getFloat(base_param+ "/density.cm3");
			}
			else
				ERRORMSG("shape for defect " << it->first << " must be one of 'disk', 'void', 'plane311' or 'irregular' and not " << aspect);
			edtype->_name = it->first;

			edtype->_transformTo = count(theArray, pPar->getString(base_param + "/to"));
			edtype->_transformFrom=count(theArray, pPar->getString(base_param + "/from"));
			edtype->_ivmodel = pPar->getBool(base_param+"/IV.model");
			edtype->_axes[0] = pPar->getCoordinates(base_param+"/axis.0");
			edtype->_axes[1] = pPar->getCoordinates(base_param+"/axis.1");
			if(edtype->_axes[0].abs() == 0 || edtype->_axes[1].abs() == 0)
				ERRORMSG("Incorrect axes in " << base_param << ". Modulus is null.");
			edtype->_doNotCreateIn = pPar->getCoordinates(base_param+"/not.in.plane");
			edtype->checkAxes();
			string migration = pPar->getString(base_param + "/migration.type");
			if(migration == "perpendicular")
				edtype->_migration = EDType::mig_perpendicular;
			else if(migration == "parallel")
				edtype->_migration = EDType::mig_parallel;
			else if(migration == "3d")
				edtype->_migration = EDType::mig_3d;
			else
				ERRORMSG("Incorrect migration in " << base_param << ". Not perpendicular, parallel or 3d");
			edtype->_lambda = pPar->getFloat(base_param+"/lambda");
			edtype->_percolation = pPar->getBool(base_param+"/percolation");
			if(pPar->specified(base_param+"/expand.impurity"))
			{
				Kernel::P_TYPE expand = pPM->getParticleNumber(mt, pPar->getString(base_param+"/expand.impurity"));
				if(expand == Kernel::UNDEFINED_TYPE)
					ERRORMSG("Incorrect expand.impurity in " << base_param);
				edtype->_expand_impurity = pPM->getFamily(expand);
				edtype->_expand_impurity_volume_nm3 = pPar->getFloat(base_param+"/expand.impurity.volume.nm3");
				edtype->_expand_capture_radius_M = pPar->getFloat(base_param+"/expand.impurity.radius.M");
				edtype->_expand_capture_radius_P = pPar->getFloat(base_param+"/expand.impurity.radius.P");
			}
			else
			{
				edtype->_expand_impurity = Kernel::UNDEFINED_TYPE;
				edtype->_expand_impurity_volume_nm3 = 0;
				edtype->_expand_capture_radius_M = 0;
				edtype->_expand_capture_radius_P = 0;
			}
			if(edtype->_transformTo >= theArray.size() || edtype->_transformFrom >= theArray.size())
				ERRORMSG("Incorrect Transform in " << base_param);
			_params[mt].push_back(edtype);  //from now we can use the cluster
			_allClusters[mt].push_back(vector<string>()); //and we can insert cluster names.
			pPar->loadProcedure(pTcl, base_param + "/formation", 0);
			pPar->loadProcedure(pTcl, base_param + "/migration", 0);
			pPar->loadProcedure(pTcl, base_param + "/prefactor", 0);
			if(edtype->_ivmodel)
				pPar->loadProcedure(pTcl, base_param + "/IV.barrier", 0);
			pPar->loadProcedure(pTcl, base_param + "/transform.to", 0);
			pPar->loadProcedure(pTcl, base_param + "/transform.from", 0);
			map<string,IO::ArrheniusAlloys> frmMap = pPar->getArrheniusAlloysProc(pTcl, base_param + "/formation", IO::Parameters::AA_ENERGY);
			map<string,IO::ArrheniusAlloys> migMap = pPar->getArrheniusAlloysProc(pTcl, base_param + "/migration", IO::Parameters::AA_FULL);
			map<string,float>               preMap = pPar->getFloatProc          (pTcl, base_param + "/prefactor");
			map<string,IO::ArrheniusAlloys> ttoMap = pPar->getArrheniusAlloysProc(pTcl, base_param + "/transform.to", IO::Parameters::AA_FULL);
			map<string,IO::ArrheniusAlloys> tfrMap = pPar->getArrheniusAlloysProc(pTcl, base_param + "/transform.from", IO::Parameters::AA_FULL);
			map<string,IO::ArrheniusAlloys> ivbMap;
			if(edtype->_ivmodel)
				ivbMap = pPar->getArrheniusAlloysProc(pTcl, base_param + "/IV.barrier", IO::Parameters::AA_FULL);
			//iterate through all of them
			for(map<string,IO::ArrheniusAlloys>::iterator it=frmMap.begin(); it!=frmMap.end(); ++it)
			{
				ID clusterID = _pPM->getID(mt, it->first);
				if(clusterID._pt.empty())
					ERRORMSG("Wrong ID in " << base_param << ": " << it->first);
				if(!isCluster(clusterID))
					ERRORMSG(base_param << ": Incorrect cluster (syntax problem or too big) " << it->first);
				if(it->first != _pPM->getIDName(clusterID, false) && it->first != _pPM->getIDName(clusterID, true))
					WARNINGMSG(base_param << ": Suggested syntax for '" << it->first << "' is '" << _pPM->getIDName(clusterID) << "'");
				_allClusters[mt].back().push_back(it->first);
				unsigned hash = _params[mt].back()->_hash.addHash(clusterID);
				//Now I have the list of clusters, so I can ask for their potentials...
				_params[mt].back()->_hash._map[hash]._eForm = it->second;
				MEDMSG(_pPM->getMaterialName(mt) << " " << " cluster " << it->first << " hash(" << hash << ") Eform=" << _params[mt].back()->_hash._map[hash]._eForm);
				//migrations
				map<string,IO::ArrheniusAlloys>::iterator itArr = migMap.find(it->first);
				if(itArr == migMap.end())
					_params[mt].back()->_hash._map[hash]._arr[0] = IO::ArrheniusAlloys(0, 5);
				else
					_params[mt].back()->_hash._map[hash]._arr[0] = itArr->second * 6./(edtype->_lambda*edtype->_lambda)*1e14;;
				//migrate correction depending on dimensionality
				if(edtype->_migration == EDType::mig_perpendicular)
					_params[mt].back()->_hash._map[hash]._arr[0] *= (2./6.);
				if(edtype->_migration == EDType::mig_parallel)
					_params[mt].back()->_hash._map[hash]._arr[0] *= (4./6.);
				MEDMSG(_pPM->getMaterialName(mt) << " " << " cluster " << it->first << " hash(" << hash <<
					") mig " << _params[mt].back()->_hash._map[hash]._arr[0]);
				//transform.to
				itArr = ttoMap.find(it->first);
				if(itArr != ttoMap.end())
					_params[mt].back()->_hash._map[hash]._arr[1] = itArr->second * convFactor;
				//transform.from
				itArr = tfrMap.find(it->first);
				if(itArr != tfrMap.end())
					_params[mt].back()->_hash._map[hash]._arr[2] = itArr->second * convFactor;
				//recombination
				itArr = ivbMap.find(it->first);
				if(itArr != ivbMap.end())
					_params[mt].back()->_hash._map[hash]._arr[3] = itArr->second * convFactor;
				else
					HIGHMSG("Barrier not found for " << it->first);
				//and emission prefactors
				for(P_TYPE pt=_pPM->getNFamilies(); pt<_pPM->getNParticles(); ++pt)
				{
					if(!_pPM->isParticleDefined(pt, mt))
						continue;
					string idx = it->first + "," + _pPM->getParticleName(mt, pt);
					map<string, float>::iterator itFlt = preMap.find(idx);
					_params[mt].back()->_hash._map[hash]._pref[pt] = (itFlt == preMap.end()? 0 : itFlt->second*convFactor);
					if(_params[mt].back()->_hash._map[hash]._pref[pt])
						MEDMSG(_pPM->getMaterialName(mt) << " cluster " << it->first << " hash(" << hash <<
						") emit " << Domains::global()->PM()->getParticleName(mt, pt) << " Pref=" << _params[mt].back()->_hash._map[hash]._pref[pt]);
				}
			}
			//syntax check of prefactors
			for(map<string,float>::iterator itCheck=preMap.begin();itCheck != preMap.end(); ++itCheck)
			{
				vector<string> tks;
				IO::ParameterManager::getTokens(itCheck->first, ',', tks);
				if(tks.size() != 2)
					ERRORMSG(base_param << ". Wrong prefactor syntax in " << itCheck->first);
				if(frmMap.find(tks[0]) == frmMap.end())
					ERRORMSG(base_param << ". Prefactor. Cannot find cluster " << tks[0] << " in formation procedure. " << itCheck->first);
				P_TYPE pt = _pPM->getParticleNumber(mt, tks[1]);
				if (pt == UNDEFINED_TYPE || !_pPM->isParticleDefined(pt, mt))
					ERRORMSG(base_param << ". Prefactor. Cannot find particle " << tks[1] << " " << itCheck->first);
			}
		}
		LOWMSG2(")");
	}
	LOWMSG("");
}

void ClusterParam::readInteractions(const IO::FileParameters * pPar)
{
	//Even the false events are to be taken, to overload possible previous true ones!
	for(M_TYPE mt=0; mt < _pPM->getNMaterials(); ++mt)
	{
		string base = _pPM->getMaterialName(mt) + "/Models/interactions";
		const float lambda = pPar->getFloat(_pPM->getMaterialName(mt) + "/Models/lambda");
		IO::array<string, string> theStrArray = pPar->getArray(base);
		for(IO::array<string, string>::iterator it=theStrArray.begin(); it!=theStrArray.end(); ++it)
		{
			vector<string> tokens_l, tokens_r;  //for the two incoming species
			IO::ParameterManager::getTokens(it->first, '+', tokens_r);
			if(tokens_r.size() != 2)
				ERRORMSG(base << " Two reactants needed for " << it->first );
			IO::ParameterManager::getTokens(tokens_r[0], ':', tokens_l);
			if(tokens_l.size() == 0 || tokens_l.size() > 2)
				ERRORMSG(base << "Unrecognized syntax in first reactant for " << it->first);

			unsigned defect_type_l = getDefectNumber(mt, tokens_l[0]);
			if(defect_type_l == _params[mt].size())
				continue;  //defect not recognize, continue;

			vector<string> tokens_result;
			IO::ParameterManager::getTokens(it->second, ',', tokens_result);
			if(tokens_result.size() != 1 && tokens_result.size() != 2)
				ERRORMSG(base << "Result of reaction not recognized in " << it->second);
			float captDist = tokens_result.size() == 2? atof(tokens_result[1].c_str()) : lambda;
			if(tokens_l.size() == 1)  // notation is ICluster+VCluster
			{
				//check that both r is a cluster
				unsigned defect_type_r = getDefectNumber(mt, tokens_r[1]);
				if(defect_type_r == _params[mt].size())
				{
					if(tokens_result[0] == "true")
					{
						Kernel::Event::E_TYPE ev = Domains::global()->PM()->getDefectType(tokens_r[1], mt);
						if(ev != Kernel::Event::INTERFACE)
							ERRORMSG(base << " Unrecognized defect type in second reactant for " << it->first);
						continue;
					}
					else
						continue;  //unrecognized cluster set to false, it is OK.
				}
				Domains::global()->setInteraction(mt, it->first);
				if(tokens_result[0] == "true")
				{
					_params[mt][defect_type_l]->_interactMC[defect_type_r] = captDist;
					_params[mt][defect_type_r]->_interactMC[defect_type_l] = captDist;
					Domains::global()->PM()->setMaximumCaptureRadius(mt, captDist);
				}
				else
				{
					_params[mt][defect_type_l]->_interactMC.erase(defect_type_r);
					_params[mt][defect_type_r]->_interactMC.erase(defect_type_l);
				}
			}
			else //notation is ICluster:I53+I
			{
				vector<ID> vecID = getIDs(mt, defect_type_l, tokens_l[1]);
				if(vecID.empty())
					WARNINGMSG(base << " Reaction " << it->first << " does not seem to math any defect");
				for(vector<ID>::const_iterator clusterID=vecID.begin(); clusterID!=vecID.end(); ++clusterID)
				{
					assert(isCluster(*clusterID)); //should be after getIDs
					unsigned hash = _params[mt][defect_type_l]->_hash.cluster2hash(*clusterID);
					if(hash == _params[mt][defect_type_l]->_hash.invalidHash())
					{
						if(tokens_result[0] != "false")
							ERRORMSG("Cluster " << tokens_l[1] << " being used without definition (" << Domains::global()->PM()->getIDName(*clusterID) << ")");
						else
							continue;  //non existing cluster set to false... it is OK
					}
					P_TYPE pt = _pPM->getParticleNumber(mt, tokens_r[1]);
					if(pt == UNDEFINED_TYPE)
					{
						if(tokens_result[0] != "false")
							ERRORMSG(base << ". Particle " << tokens_r[1] << " could not be recognized!");
						else
							continue;
					}
					else
					{
						if(!_pPM->isParticleDefined(pt, mt))
							continue;
						if(tokens_result[0] == "sink")
						{
							_params[mt][defect_type_l]->_hash._map[hash]._interactMP_sink[pt] = true;
							_params[mt][defect_type_l]->_hash._map[hash]._interactMP[pt] = true;
							_params[mt][defect_type_l]->_hash._map[hash]._interactMP_radius[pt] = captDist;
							Domains::global()->PM()->setMaximumCaptureRadius(mt, captDist);
						}
						else
						{
							//check that result exists
							Kernel::ID ID2 = pt2ID(mt, pt);
							Kernel::ID ID1 = *clusterID;
							Cluster::addMap(ID1, ID2, _params[mt][defect_type_l]->_ivmodel);
							unsigned newHash = _params[mt][defect_type_l]->_hash.cluster2hash(ID1);
							P_TYPE ppt = ID2pt(ID1);
							if( newHash != _params[mt][defect_type_l]->_hash.invalidHash() ||
								(ppt != UNDEFINED_TYPE && _pPM->isParticleDefined(ppt, mt))	)
							{
								if(tokens_result[0] != "false")
								{
									_params[mt][defect_type_l]->_hash._map[hash]._interactMP_sink[pt] = false;
									_params[mt][defect_type_l]->_hash._map[hash]._interactMP[pt] = true;
									_params[mt][defect_type_l]->_hash._map[hash]._interactMP_radius[pt] = captDist;
									Domains::global()->PM()->setMaximumCaptureRadius(mt, captDist);
								}
								else
								{
									_params[mt][defect_type_l]->_hash._map[hash]._interactMP[pt] = false;
									_params[mt][defect_type_l]->_hash._map[hash]._interactMP_sink[pt] = false;
								}
								MEDMSG(_pPM->getMaterialName(mt) << " " << tokens_r[0] << "+" << _pPM->getParticleName(mt, pt));
								Domains::global()->setInteraction(mt, it->first);
							}
						}
					}
				} //end for vector<ID>
			}
		}
	}
}

//a cluster: Is not a particle, and does not contain at the same type Is and Vs.
bool ClusterParam::isCluster(const ID &theMap)
{
	if(ID2pt(theMap) != UNDEFINED_TYPE)
		return false;
	return true;
}

P_TYPE ClusterParam::ID2pt(const ID &theMap)
{
	P_TYPE mt0 = Domains::global()->PM()->getMaterial(theMap._mt)._pt[0];
	P_TYPE mt1 = Domains::global()->PM()->getMaterial(theMap._mt)._pt[1];
	const bool binary = Domains::global()->PM()->getMaterial(theMap._mt)._binary;
	unsigned nDop = 0;
	P_TYPE theDop = UNDEFINED_TYPE;
	unsigned nPos[NO_POS] = { 0, 0, 0, 0 };
	unsigned nV = 0, nPt1 = 0, nPt2 = 0;
	for(map<P_TYPE, unsigned>::const_iterator it=theMap._pt.begin(); it!=theMap._pt.end(); ++it)
	{
		if(it->first >= Domains::global()->PM()->getNFamilies())
			ERRORMSG(Domains::global()->PM()->getRawIDName(theMap) << " fam: Incorrect ID detected in ClusterParam::ID2pt for " << Domains::global()->PM()->getParticleName(theMap._mt, it->first));
		if(it->first == V_TYPE)
			nV+= it->second;
		else if(it->first == mt0)
			nPt1 += it->second;
		else if(it->first == mt1)
			nPt2 += it->second;
		else
		{
			theDop = it->first;
			nDop+= it->second;
		}
	}
	for(map<P_POS, unsigned>::const_iterator it=theMap._pos.begin(); it!=theMap._pos.end(); ++it)
		nPos[it->first] += it->second;
	if(nPos[POS_V])
		ERRORMSG("ID2pt: Incorrect ID at pos.");
	unsigned total = nDop + nV + nPt1 + nPt2;
	if(total == 1)
	{
		P_TYPE pt = Domains::global()->PM()->getParticle(theMap._mt, theMap._pt.begin()->first, theMap._pos.begin()->first);
		if(pt == UNDEFINED_TYPE)
			return pt;
		return (Domains::global()->PM()->isParticleDefined(pt, theMap._mt) ? pt : UNDEFINED_TYPE);
	}
	if(total == 2 && nV == 1 && nDop == 1) //HeV, etc...
	{
		if((binary && (nPos[POS_0] != 1 || nPos[POS_1] != 1)) ||
		   (!binary && nPos[POS_0] != 2))
			ERRORMSG("Wrong ID in Cluster::ID2pt");
		P_TYPE pt = Domains::global()->PM()->getParticle(theMap._mt, theDop, POS_V);
		return (Domains::global()->PM()->isParticleDefined(pt, theMap._mt) ? pt : UNDEFINED_TYPE);
	}
	return UNDEFINED_TYPE;
}


const EDType::CLType * ClusterParam::getParams(unsigned et, const ID &m) const
{
	if(m._pt.size() == 0)
		return 0;
	unsigned hash = _params[m._mt][et]->_hash.cluster2hash(m);
	return getParams(m._mt, et, hash);
}

const EDType::CLType * ClusterParam::getParams(Kernel::M_TYPE mt, unsigned et, unsigned hash) const
{
	if(hash == _params[mt][et]->_hash.invalidHash())
		return 0;
	return &_params[mt][et]->_hash._map[hash];
}

//understands templates like B*I*
vector<Kernel::ID> ClusterParam::getIDs(M_TYPE mt, unsigned edtype, const std::string &txt) const
{
	vector<Kernel::ID> storage;
	for(vector<std::string>::const_iterator it=_allClusters[mt][edtype].begin(); it!=_allClusters[mt][edtype].end(); ++it)
		if(fnmatch(txt.c_str(), it->c_str(), 0) == 0) //if matches
			storage.push_back(_pPM->getID(mt, *it));
	return storage;
}


unsigned ClusterParam::count(const map<string,bool> &array, const string &key)
{
	unsigned c=0;
	for(map<string, bool>::const_iterator it=array.begin(); it!=array.end(); ++it)
	{
		if(!it->second)
			continue;
		if(it->first == key)
			return c;
		c++;
	}
	return c;
}

unsigned ClusterParam::getDefectNumber(Kernel::M_TYPE mt, const std::string &key) const
{
	for(unsigned i=0; i<_params[mt].size(); ++i)
			if(_params[mt][i]->_name == key)
				return i;
	return _params[mt].size();
}


bool ClusterParam::reactionPossible(Kernel::SubDomain *pSub, Kernel::Domain *pDom, const Kernel::MeshElement *pEle, unsigned edtype, P_TYPE pt1, P_TYPE pt2, float kT) const
{
	Kernel::M_TYPE mt = pEle->getMaterial();
	ID theMap = pt2ID(mt, pt1);
	Cluster::addTo(pt2, theMap, _params[mt][edtype]->_ivmodel);
	const EDType::CLType *otherCl = getParams(edtype, theMap);
	if(otherCl == 0)
		return false; //no reaction possible
	float init = pDom->_pMPPar->_form[mt][pt1](pEle)._ener + pDom->_pMPPar->_form[mt][pt2](pEle)._ener;
	float end  = otherCl->_eForm(pEle)._ener;
	double arg = (init - end);
	return (arg < 0? pSub->_rng.rand() < exp(arg /kT) : true);
}

ID ClusterParam::pt2ID(Kernel::M_TYPE mt, P_TYPE pt)
{
	P_TYPE fam = Domains::global()->PM()->getFamily(pt);
	P_POS  pos = Domains::global()->PM()->getPPos(pt);
	bool binary = Domains::global()->PM()->getMaterial(mt)._binary;
	ID id;
	id._mt = mt;
	if(pos == POS_V) //a pair, actually
	{
		if(binary)
		{
			id._pos[POS_0] = 1;
			id._pos[POS_1] = 1;
		}
		else
			id._pos[POS_0] = 2;
		id._pt[V_TYPE] = 1;
	}
	else
		id._pos[pos] = 1;
	id._pt[fam] = 1;
	return id;
}

}

